#!/usr/bin/env python3


"""Classify splice sites with a neural network and SVM.

Usage:
  classify-splice-sites.py <data-file-name> <model-file-name> <window-inner> \
    <window-outer> <site> [-r <imbalance-ratio>] [-c <ncpus>]
  classify-splice-sites.py (-h | --help)

Options:
  -h --help     Show this screen.
  -c <ncpus>    Number of CPUs  [default: 4]
"""

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from docopt import docopt

from tools import performance_metrics
from tools import read_data
from tools import read_model

DNA_SYMBOLS = np.array(['A', 'T', 'C', 'G', 'N'])


def main():
    arguments = docopt(__doc__, version='1.0')

    data_file = Path(arguments['<data-file-name>'])
    model_file = arguments['<model-file-name>']

    window_inner = int(arguments['<window-inner>'])
    window_outer = int(arguments['<window-outer>'])

    site = arguments['<site>']
    assert site in ('donor', 'acceptor')

    ncpus = int(arguments['-c'])

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        filename=f'classify-splice-sites-{os.path.basename(model_file)}.log',
        filemode='w'
    )
    logging.info(f'Classification of file {data_file} with model {model_file}')

    model_extension = os.path.splitext(model_file)[1]
    if model_extension == '.hd5':
        input_df, predictions = get_svm_predictions(model_file, data_file, window_inner, window_outer, site, ncpus)
    elif model_extension == '.h5':
        input_df, predictions = get_dnn_predictions(model_file, data_file, window_inner, window_outer, site)
    else:
        logging.warning(f'Unknown model extension {model_extension}. Exiting.')
        return

    output_df = input_df.assign(pred=pd.Series(predictions))
    output_df.to_csv(sys.stdout, sep=';', index=False)

    if 'label' in input_df:
        imbalance_ratio = float(arguments['<imbalance-ratio>'])

        log_performance(input_df[["label"]], predictions, imbalance_ratio)


def prepare_inputs(input_df: pd.DataFrame, site: str):
    sequences = []
    for sequence in input_df['sequence']:
        sequence = str(sequence)
        sequence = sequence.upper()

        try:
            assert set(sequence).issubset({'A', 'C', 'T', 'G', 'N'})
        except AssertionError:
            print(set(sequence))

        if site == 'acceptor':
            assert sequence[198] == 'A'
            assert sequence[199] == 'G'
        else:
            assert sequence[200] == 'G'
            assert sequence[201] == 'T'

        sequence = np.array(list(sequence))

        sequences.append((sequence[:, None] == DNA_SYMBOLS).astype(np.float32))

    return np.array(sequences)


def get_svm_predictions(model_file, data_file, window_inner, window_outer, site, ncpus):
    window = (window_outer, window_inner) if site == 'donor' else (window_inner, window_outer)

    # Extract the region of interest given the window.
    # NOTE: The script expects the acceptor/donor dimer to be in the middle of the sequence
    input_df = read_data(data_file, window=window)

    import shogun as sg  # import it here not to conflict with keras d-NN use
    sg.Parallel().set_num_threads(ncpus)

    features = sg.StringCharFeatures(input_df.sequence.tolist(), sg.RAWBYTE)
    model = read_model(model_file)

    predict = model.apply_binary(features)
    return input_df, predict.get_int_labels()


def get_dnn_predictions(model_file, data_file, window_inner, window_outer, site):
    assert window_inner == 200
    assert window_outer == 200

    window = (window_outer, window_inner - 2) if site == 'donor' else (window_inner - 2, window_outer)

    os.environ['KERAS_BACKEND'] = 'tensorflow'
    from tensorflow.keras.models import load_model

    input_df = read_data(data_file, window=window)
    model = load_model(model_file)

    inputs = prepare_inputs(input_df, site)
    predictions = model.predict(inputs)
    predictions = np.squeeze(predictions)
    predictions[np.where(predictions < 0.65)] = -1  # negative classes are -1 as opposed to NN output (which is 0)
    predictions[np.where(predictions >= 0.65)] = 1
    predictions = predictions.astype(np.int32)

    return input_df, predictions


def log_performance(prediction, labels, imbalance_ratio):
    if imbalance_ratio is None:
        # Set a new reasonable default (but still wrong) value
        imbalance_ratio = 0.005
        logging.warning(f'Imbalance ratio was not set! Using default value {imbalance_ratio}. '
                        f'Results will be inaccurate!')

    metrics_list = performance_metrics(prediction, labels, imbalance_ratio)
    metrics_str = '\n'.join(metrics_list)
    logging.info(f'Performance metrics:\n'
                 f'{metrics_str}')


if __name__ == '__main__':
    main()
