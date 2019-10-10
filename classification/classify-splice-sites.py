import sys
from contextlib import closing

import pandas as pd
import shogun as sg


def read_model(filename):
    svm = sg.LibSVM()

    model_file = sg.SerializableHdf5File(filename, "r")

    with closing(model_file):
        if not svm.load_serializable(model_file):
            exit(1)

    return svm


def parser():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Classify splice sites")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('model_filename', metavar='MODEL', type=str,
                   help='filename of the model')
    p.add_argument('-c', '--cpus', type=int, default=1,
                   dest='ncpus',
                   help='number of CPUs')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    sg.Parallel().set_num_threads(argparser.ncpus)

    data = pd.read_csv(argparser.data_filename, sep=';')
    model = read_model(argparser.model_filename)

    features = sg.StringCharFeatures(data.sequence.tolist(), sg.RAWBYTE)

    predict = model.apply_binary(features)

    data.assign(pred=pd.Series(list(predict.get_int_labels()))) \
        .to_csv(sys.stdout, sep=';', index=False)
