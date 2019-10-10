from contextlib import closing

import numpy as np
import pandas as pd
import shogun as sg
from sklearn.model_selection import StratifiedShuffleSplit


def split_data(data, test_size):
    s = StratifiedShuffleSplit(n_splits=1, test_size=test_size)
    for train_i, test_i in s.split(data.sequence, data.label):
        train, test = data.loc[train_i], data.loc[test_i]
        return train, test


def train(sequences, labels, d, C_n, C_p=None):
    if C_p is None:
        print("Training SVM(C={}, d={})".format(C_n, d))
    else:
        print("Training SVM(C_neg={}, C_pos={}, d={})".format(C_n, C_p, d))

    print()

    features = sg.StringCharFeatures(sequences.tolist(), sg.RAWBYTE)
    labels = sg.BinaryLabels(np.array(labels))

    kernel = sg.WeightedDegreeStringKernel(features, features, d)
    svm = sg.LibSVM(C_n, kernel, labels)

    if C_p is not None:
        svm.set_C(C_n, C_p)

    svm.train()

    predict = svm.apply_binary(features)

    acc = sg.AccuracyMeasure()
    acc.evaluate(predict, labels)
    print("Train results:")
    print("TP", "FP", "TN", "FN", sep='\t')
    print(int(acc.get_TP()), int(acc.get_FP()),
          int(acc.get_TN()), int(acc.get_FN()), sep='\t')

    print()

    return svm


def test(model, sequences, labels):
    features = sg.StringCharFeatures(sequences.tolist(), sg.RAWBYTE)
    labels = sg.BinaryLabels(np.array(labels))

    predict = model.apply_binary(features)

    acc = sg.AccuracyMeasure()
    acc.evaluate(predict, labels)

    print("Test results:")
    print("TP", "FP", "TN", "FN", sep='\t')
    print(int(acc.get_TP()), int(acc.get_FP()),
          int(acc.get_TN()), int(acc.get_FN()), sep='\t')

    print()


def parser():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Train splice site model")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('d', metavar='DEGREE', type=int,
                   help='degree of the WD kernel')
    p.add_argument('C', metavar='C', type=float,
                   nargs='+',
                   help='the regularization constant, if two values are given then it is considered as C_neg C_pos')
    p.add_argument('-o', '--output', type=str, default="model.hd5",
                   dest='model_filename',
                   help='filename for the trained model (default is \'model.hd5\')')
    p.add_argument('-t', '--test', type=float, default=0,
                   dest='test_size',
                   help='what fraction of the data should be used for testing, if not specified then only training will be performed')
    p.add_argument('-c', '--cpus', type=int, default=1,
                   dest='ncpus',
                   help='number of CPUs')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    sg.Parallel().set_num_threads(argparser.ncpus)

    data = pd.read_csv(argparser.data_filename, sep=';')

    if argparser.test_size > 0:
        train_data, test_data = split_data(data, argparser.test_size)
    else:
        train_data = data

    model = train(train_data.sequence, train_data.label,
                  argparser.d, *argparser.C)

    if argparser.test_size > 0:
        test(model, test_data.sequence, test_data.label)

    model_file = sg.SerializableHdf5File(argparser.model_filename, "w")
    with closing(model_file):
        if model.save_serializable(model_file):
            print("Model saved: {}".format(argparser.model_filename))
        else:
            print("Model could not be saved")
