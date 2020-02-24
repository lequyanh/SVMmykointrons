#!/usr/bin/env python3


"""Classify introns with SVM.

Usage:
  classify-introns.py <data-file-name> <model-file-name> <kernel-order> [-c <ncpus>]
  classify-introns.py (-h | --help)

Options:
  -h --help     Show this screen.
  -c <ncpus>    Number of CPUs  [default: 4]
"""

import logging
import sys

import numpy as np
import pandas as pd
import shogun as sg
from docopt import docopt

from tools import performance_metrics, read_model


def create_features(order, dna):
    gap = 0
    reverse = False

    charfeat = sg.StringCharFeatures(sg.DNA)
    charfeat.set_features(dna)
    feats = sg.StringWordFeatures(charfeat.get_alphabet())
    feats.obtain_from_char(charfeat, order - 1, order, gap, reverse)
    preproc = sg.SortWordString()
    preproc.init(feats)
    return preproc.apply(feats)


if __name__ == "__main__":
    arguments = docopt(__doc__, version='1.0')

    data_file = arguments['<data-file-name>']
    model_file = arguments['<model-file-name>']

    kernel_order = int(arguments['<kernel-order>'])  # order of the spectrum kernel
    ncpus = int(arguments['-c'])

    data = pd.read_csv(data_file, sep=';')
    model = read_model(model_file)

    sg.Parallel().set_num_threads(ncpus)

    features = create_features(kernel_order, data['sequence'].tolist())
    predict = model.apply_binary(features)

    data.assign(prediction=pd.Series(list(predict.get_int_labels()))) \
        .to_csv(sys.stdout, sep=';', index=False)

    if 'label' in data:
        logging.basicConfig(
            level=logging.INFO,
            filename=f'classify-introns-d{kernel_order}.log',
            filemode='w'
        )
        logging.info(f'Validation of model {model_file} on testset {data_file}; d={kernel_order}')

        labels = sg.BinaryLabels(np.array(data.label))

        metrics_data = performance_metrics(
            labels.get_int_labels(),
            predict.get_int_labels(),
            imbalance_rat=None
        )

        metrics_str = '\n'.join(metrics_data)
        logging.info(f'Performance metrics:\n'
                     f'{metrics_str}')
