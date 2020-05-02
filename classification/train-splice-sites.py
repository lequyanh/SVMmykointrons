#!/usr/bin/env python3


"""Classify splice sites with a neural network and SVM.

Usage:
  train-splice-sites.py <data-file-name> <model-file-name> <window-inner> \
    <window-outer> <site> <kernel-degree> <regularization-const> \
    [-t <test-size>] [-r <imbalance-ratio>] [-c <ncpus>]
  train-splice-sites.py (-h | --help)

Options:
  -h --help         Show this screen.
  -c <ncpus>        Number of CPUs  [default: 4]
  -t <test-size>    Test size  [default: 0]
"""

import logging
from contextlib import closing

import numpy as np
import shogun as sg
from docopt import docopt

from tools import read_data, split_data, performance_metrics


def main():
    arguments = docopt(__doc__, version='1.0')

    data_file = arguments['<data-file-name>']
    out_model_name = arguments['<model-file-name>']

    window_inner = int(arguments['<window-inner>'])
    window_outer = int(arguments['<window-outer>'])

    site = arguments['<site>']
    assert site in ('donor', 'acceptor')

    kernel_degree = int(arguments['<kernel-degree>'])
    C = int(arguments['<regularization-const>'])

    test_size = float(arguments['-t'])
    imbalance_ratio = float(arguments['<imbalance-ratio>'])

    number_of_cpus = int(arguments['-c'])

    window = (window_outer, window_inner) if site == 'donor' else (window_inner, window_outer)
    data = read_data(data_file, window=window)

    logging.basicConfig(
        level=logging.INFO,
        filename=f'train-splice-sites-d{kernel_degree}-C{C}-{window_inner}-{window_outer}.log',
        filemode='w'
    )
    logging.info(f'Model of trained on {data_file}')

    if imbalance_ratio is None:
        # Set a new reasonable default (but still wrong) value
        imbalance_ratio = 0.005
        logging.warning(f'Imbalance ratio was not set! Using default value {imbalance_ratio}. '
                        f'Adjusted precision will be incorrect!')

    if test_size > 0:
        train_data, test_data = split_data(data, test_size)

        train_features = sg.StringCharFeatures(train_data.sequence.tolist(), sg.RAWBYTE)
        test_features = sg.StringCharFeatures(test_data.sequence.tolist(), sg.RAWBYTE)

        train_labels = sg.BinaryLabels(np.array(train_data.label))
        test_labels = np.array(test_data.label)

    else:
        train_data = data

        train_features = sg.StringCharFeatures(train_data.sequence.tolist(), sg.RAWBYTE)
        train_labels = sg.BinaryLabels(np.array(train_data.label))

    sg.Parallel().set_num_threads(number_of_cpus)
    svm_model = train(train_features, train_labels, kernel_degree, C, imbalance_ratio)

    if test_size > 0:
        test(svm_model, test_features, test_labels, imbalance_ratio)

    output = sg.SerializableHdf5File(out_model_name, "w")
    with closing(output):
        if svm_model.save_serializable(output):
            print("Model saved: {}".format(out_model_name))
        else:
            print("Saving the model failed")


def train(train_features, train_labels, kernel_degree, C, imbalance_ratio):
    wd_kernel = sg.WeightedDegreeStringKernel(train_features, train_features, kernel_degree)

    svm_model = sg.LibSVM(C, wd_kernel, train_labels)
    svm_model.train()

    prediction = svm_model.apply_binary(train_features)

    metrics_list = performance_metrics(train_labels.get_int_labels(), prediction.get_int_labels(), imbalance_ratio)
    metrics_str = '\n'.join(metrics_list)
    logging.info(f'Train performance metrics:\n'
                 f'{metrics_str}\n')

    return svm_model


def test(svm_model, test_features, test_labels, imbalance_ratio):
    prediction = svm_model.apply_binary(test_features)
    prediction = prediction.get_int_labels()

    metrics_list = performance_metrics(test_labels, prediction, imbalance_ratio)
    metrics_str = '\n'.join(metrics_list)
    logging.info(f'Test performance metrics:\n'
                 f'{metrics_str}')


if __name__ == "__main__":
    main()
