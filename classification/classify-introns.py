import logging
import sys
from contextlib import closing

import numpy as np
import pandas as pd
import shogun as sg


def read_model(filename):
    svm = sg.LibSVM()

    model_file = sg.SerializableHdf5File(filename, "r")

    with closing(model_file):
        if not svm.load_serializable(model_file):
            print("Model failed to load")
            exit(1)

    return svm


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


def parser():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Classify introns")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('model_filename', metavar='MODEL', type=str,
                   help='filename of the model')
    p.add_argument('l', metavar='ORDER', type=int,
                   help='order of the spectrum kernel')
    p.add_argument('-c', '--cpus', type=int, default=1,
                   dest='ncpus',
                   help='number of CPUs')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    sg.Parallel().set_num_threads(argparser.ncpus)

    data = pd.read_csv(argparser.data_filename, sep=';')
    model = read_model(argparser.model_filename)

    features = create_features(argparser.l, data.loc[:, 'sequence'].tolist())

    predict = model.apply_binary(features)

    data.assign(pred=pd.Series(list(predict.get_int_labels()))) \
        .to_csv(sys.stdout, sep=';', index=False)

    if 'label' in data:
        logging.basicConfig(
            level=logging.INFO,
            filename=f'classification-d{argparser.l}.log',
            filemode='w'
        )

        logging.info(
            f'Validation of model {argparser.model_filename} on testset {argparser.data_filename}; d={argparser.l}'
        )

        labels = sg.BinaryLabels(np.array(data.label))
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
