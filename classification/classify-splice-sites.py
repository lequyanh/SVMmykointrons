import logging
import os
import sys
from tools import read_data, read_model

import numpy as np
import pandas as pd
import shogun as sg


def parser():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Classify splice sites")
    p.add_argument('data_filename', metavar='INPUT', type=str,
                   help='filename of the input')
    p.add_argument('model_filename', metavar='MODEL', type=str,
                   help='filename of the model')
    p.add_argument('window_inner', type=int)
    p.add_argument('window_outer', type=int)
    p.add_argument('site', type=str, help='donor or acceptor')
    p.add_argument('-v', action='store_true', help='validation mode - calculates accuracy metrics')
    p.add_argument('-c', '--cpus', type=int, default=1,
                   dest='ncpus',
                   help='number of CPUs')
    return p


if __name__ == "__main__":
    argparser = parser().parse_args()

    sg.Parallel().set_num_threads(argparser.ncpus)

    if argparser.v:
        win_in = argparser.window_inner
        win_out = argparser.window_outer
        window = (win_out, win_in) if argparser.site == 'donor' else (win_in, win_out)

        data = read_data(argparser.data_filename, window=window)
    else:
        data = pd.read_csv(argparser.data_filename, sep=';')

    features = sg.StringCharFeatures(data.sequence.tolist(), sg.RAWBYTE)
    model = read_model(argparser.model_filename)

    predict = model.apply_binary(features)

    if argparser.v:
        mod_name = os.path.basename(argparser.model_filename)
        logging.basicConfig(
            level=logging.INFO,
            filename=f'classify-splice-sites-{mod_name}.log',
            filemode='w'
        )

        logging.info(f'Classification of file {argparser.data_filename} with model {argparser.model_filename}')

        labels = sg.BinaryLabels(np.array(data.label))
        acc = sg.AccuracyMeasure()
        metrics = ["Accuracy: {}".format(acc.evaluate(predict, labels)),
                   " -TP: {}".format(acc.get_TP()),
                   " -FP: {}".format(acc.get_FP()),
                   " -TN: {}".format(acc.get_TN()),
                   " -FN: {}".format(acc.get_FN())]
        metrics = '\n'.join(metrics)
        logging.info(metrics)

    data.assign(pred=pd.Series(list(predict.get_int_labels()))) \
        .to_csv(sys.stdout, sep=';', index=False)
