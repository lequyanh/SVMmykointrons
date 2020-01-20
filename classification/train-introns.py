import logging
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


def create_features(order, train_dna, test_dna=None):
    gap = 0
    reverse = False

    charfeat = sg.StringCharFeatures(sg.DNA)
    charfeat.set_features(train_dna)
    feats_train = sg.StringWordFeatures(charfeat.get_alphabet())
    feats_train.obtain_from_char(charfeat, order - 1, order, gap, reverse)
    preproc = sg.SortWordString()
    preproc.init(feats_train)
    feats_train = preproc.apply(feats_train)

    if test_dna:
        charfeat = sg.StringCharFeatures(sg.DNA)
        charfeat.set_features(test_dna)
        feats_test = sg.StringWordFeatures(charfeat.get_alphabet())
        feats_test.obtain_from_char(charfeat, order - 1, order, gap, reverse)
        feats_test = preproc.apply(feats_test)

        return feats_train, feats_test
    else:
        return feats_train


def train(features, labels, cache, l, C_n, C_p=None):
    if C_p is None:
        logging.info("Training SVM(C={}, l={})".format(C_n, l))
    else:
        logging.info("Training SVM(C_neg={}, C_pos={}, l={})".format(C_n, C_p, l))

    kernel = sg.CommWordStringKernel(features, features, False, cache)
    svm = sg.LibSVM(C_n, kernel, labels)

    if C_p is not None: svm.set_C(C_n, C_p)

    svm.train()

    predict = svm.apply_binary(features)

    acc = sg.AccuracyMeasure()
    acc.evaluate(predict, labels)
    TP = int(acc.get_TP())
    FP = int(acc.get_FP())
    FN = int(acc.get_FN())
    TN = int(acc.get_TN())

    logging.info("Train results:")
    logging.info('\t'.join(["TP", "FP", "TN", "FN"]))
    logging.info('\t'.join(map(str, [TP, FP, TN, FN])))

    logging.info(f'Precision: {TP / (TP + FP)}\n'
                 f'Recall: {TP / (TP + FN)}\n'
                 f'Accuracy: {(TP + TN) / (sum([TP, FP, FN, TN]))}')

    logging.info("")

    return svm


def test(model, features, labels):
    predict = model.apply_binary(features)

    acc = sg.AccuracyMeasure()
    acc.evaluate(predict, labels)
    TP = int(acc.get_TP())
    FP = int(acc.get_FP())
    FN = int(acc.get_FN())
    TN = int(acc.get_TN())

    logging.info("Test results:")
    logging.info('\t'.join(["TP", "FP", "TN", "FN"]))
    logging.info('\t'.join(map(str, [TP, FP, TN, FN])))

    logging.info(f'Precision: {TP / (TP + FP)}\n'
                 f'Recall: {TP / (TP + FN)}\n'
                 f'Accuracy: {(TP + TN) / (sum([TP, FP, FN, TN]))}')


def parser():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Train intron model")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('l', metavar='ORDER', type=int,
                   help='order of the spectrum kernel')
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
    p.add_argument('-s', '--cache', type=int, default=1024,
                   dest='cache_size',
                   help='size of cache used by the spectrum kernel')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    logging.basicConfig(
        level=logging.INFO,
        filename=f'train-introns-d{argparser.l}-C{argparser.C}.output',
        filemode='w'
    )

    sg.Parallel().set_num_threads(argparser.ncpus)

    data = pd.read_csv(argparser.data_filename, sep=';')

    if argparser.test_size > 0:
        train_data, test_data = split_data(data, argparser.test_size)
        train_features, test_features = create_features(argparser.l,
                                                        train_data.loc[:, 'sequence'].tolist(),
                                                        test_data.loc[:, 'sequence'].tolist())
    else:
        train_data = data
        train_features = create_features(argparser.l, train_data.loc[:, 'sequence'].tolist())

    model = train(train_features, sg.BinaryLabels(np.array(train_data.label)),
                  argparser.cache_size, argparser.l, *argparser.C)

    if argparser.test_size > 0:
        test(model, test_features, sg.BinaryLabels(np.array(test_data.label)))

    model_file = sg.SerializableHdf5File(argparser.model_filename, "w")
    with closing(model_file):
        if model.save_serializable(model_file):
            logging.info("Model saved: {}".format(argparser.model_filename))
        else:
            logging.warning("Model could not be saved")
