import logging
import os
import sys
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import shogun as sg

from tools import read_data, read_model, performance_metrics


def parser():
    p = ArgumentParser(description="Classify splice sites")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('model_filename', metavar='MODEL', type=str,
                   help='filename of the model')
    p.add_argument('window_inner', type=int)
    p.add_argument('window_outer', type=int)
    p.add_argument('site', type=str, help='donor or acceptor')
    p.add_argument('-r', type=float, default=-1.0,
                   help='real imbalance ratio of +/- classes in data. Float in <0,1>')
    p.add_argument('-c', '--cpus', type=int, default=1, dest='ncpus')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    win_in = argparser.window_inner
    win_out = argparser.window_outer
    window = (win_out, win_in) if argparser.site == 'donor' else (win_in, win_out)

    # Extract the region of interest given the window.
    # NOTE: The script expects the acceptor/donor dimer to be in the middle of the sequence
    data = read_data(argparser.data_filename, window=window)

    sg.Parallel().set_num_threads(argparser.ncpus)

    features = sg.StringCharFeatures(data.sequence.tolist(), sg.RAWBYTE)
    model = read_model(argparser.model_filename)

    predict = model.apply_binary(features)

    if 'label' in data:
        mod_name = os.path.basename(argparser.model_filename)
        logging.basicConfig(
            level=logging.INFO,
            filename=f'classify-splice-sites-{mod_name}.log',
            filemode='w'
        )
        logging.info(f'Classification of file {argparser.data_filename} with model {argparser.model_filename}')

        labels = sg.BinaryLabels(np.array(data.label))
        imbalance_ratio = argparser.r

        if imbalance_ratio < 0:
            # Dummy value was set by argparser, set a new reasonable default (but still wrong) value
            imbalance_ratio = 0.005
            logging.warning(f'Imbalance ratio was not set! Using default value {imbalance_ratio}. '
                            f'Results will be inaccurate!')

        metrics_list = performance_metrics(predict, labels, imbalance_ratio)
        metrics_str = '\n'.join(metrics_list)
        logging.info(f'Performance metrics:\n'
                     f'{metrics_str}')

    data.assign(
        pred=pd.Series(list(predict.get_int_labels()))) \
        .to_csv(sys.stdout, sep=';', index=False)
