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
import shogun as sg
from docopt import docopt

from tools import performance_metrics
from tools import read_data
from tools import read_model

DNA_SYMBOLS = np.array(['A', 'T', 'C', 'G', 'N'])
BASIDIO_DONOR_FREQ = 1 / 40
np.random.seed(42)


def main():
    arguments = docopt(__doc__, version='1.0')

    data_file = Path(arguments['<data-file-name>'])
    model_file = arguments['<model-file-name>']

    window_inner = int(arguments['<window-inner>'])
    window_outer = int(arguments['<window-outer>'])

    site = arguments['<site>']
    assert site in ('donor', 'acceptor')

    ncpus = int(arguments['-c'])

    if arguments['<imbalance-ratio>']:
        imbalance_ratio = float(arguments['<imbalance-ratio>'])
        print(f'Imbalance ratio {imbalance_ratio}')

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        filename=f'classify-splice-sites-{os.path.basename(model_file)}.log',
        filemode='w'
    )
    logging.info(f'Classification of file {data_file} with model {model_file}')

    if model_file == 'random':
        input_df, predictions = get_random_predictions(data_file, window_inner, window_outer, site)
    else:
        model_extension = os.path.splitext(model_file)[1]
        if model_extension == '.hd5':
            input_df, predictions = get_svm_predictions(model_file, data_file, window_inner, window_outer, site, ncpus)
        else:
            logging.warning(f'Unknown model extension {model_extension}. Exiting.')
            return

    output_df = input_df.assign(pred=pd.Series(predictions))
    output_df.to_csv(sys.stdout, sep=';', index=False)

    if 'label' in input_df:
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

        if site == 'acceptor' and len(sequence) == 400:
            assert sequence[198] == 'A'
            assert sequence[199] == 'G'
        elif site == 'donor' and len(sequence) == 400:
            assert sequence[200] == 'G'
            assert sequence[201] == 'T'
        elif site == 'acceptor' and len(sequence) == 100:
            assert sequence[98] == 'A'
            assert sequence[99] == 'G'
        elif site == 'donor' and len(sequence) == 100:
            assert sequence[0] == 'G'
            assert sequence[1] == 'T'

        sequence = np.array(list(sequence))

        sequences.append((sequence[:, None] == DNA_SYMBOLS).astype(np.float32))

    return np.array(sequences)


def get_svm_predictions(model_file, data_file, window_inner, window_outer, site, ncpus):
    window = (window_outer, window_inner) if site == 'donor' else (window_inner, window_outer)

    # Extract the region of interest given the window.
    # NOTE: The script expects the acceptor/donor dimer to be in the middle of the sequence
    input_df = read_data(data_file, window=window)

    logging.info(f'{sg.Parallel().get_num_threads()} threads')

    features = sg.StringCharFeatures(input_df.sequence.tolist(), sg.RAWBYTE)
    model = read_model(model_file)

    predict = model.apply_binary(features)
    return input_df, predict.get_int_labels()


def get_random_predictions(data_file, window_inner, window_outer, site, ):
    window = (window_outer, window_inner) if site == 'donor' else (window_inner, window_outer)

    input_df = read_data(data_file, window)  # window can be random
    n = len(input_df)
    f = 2 * BASIDIO_DONOR_FREQ

    frequencies = [1 - f, f]

    predictions = np.random.choice([0, 1], n, replace=True, p=frequencies)

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
