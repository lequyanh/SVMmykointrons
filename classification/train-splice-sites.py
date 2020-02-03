import logging
import sys
from contextlib import closing

import numpy as np
import shogun as sg

from tools import read_data

if __name__ == "__main__":
    data_file = sys.argv[1]
    model_file = sys.argv[2]

    lWindow = int(sys.argv[3])
    rWindow = int(sys.argv[4])

    kernel_degree = int(sys.argv[5])
    C = int(sys.argv[6])

    number_of_cpus = int(sys.argv[7])

    sg.Parallel().set_num_threads(number_of_cpus)

    data = read_data(data_file, window=(lWindow, rWindow))

    features = sg.StringCharFeatures(data.sequence.tolist(), sg.RAWBYTE)
    labels = sg.BinaryLabels(np.array(data.label))

    wd_kernel = sg.WeightedDegreeStringKernel(features, features, kernel_degree)

    svm = sg.LibSVM(C, wd_kernel, labels)
    svm.train()

    predict = svm.apply_binary(features)

    acc = sg.AccuracyMeasure()

    ACC = acc.evaluate(predict, labels)
    TP = int(acc.get_TP())
    FP = int(acc.get_FP())
    FN = int(acc.get_FN())
    TN = int(acc.get_TN())

    logging.basicConfig(
        level=logging.INFO,
        filename=f'train-splice-sites-d{kernel_degree}-C{C}-{lWindow}-{rWindow}.log',
        filemode='w'
    )

    logging.info(f'Model of trained on {data_file}')

    logging.info("Train results:")
    logging.info('\t'.join(["TP", "FP", "TN", "FN"]))
    logging.info('\t'.join(map(str, [TP, FP, TN, FN])))

    logging.info(f'Precision: {TP / (TP + FP)}\n'
                 f'Recall: {TP / (TP + FN)}\n'
                 f'Accuracy: {ACC}')

    output = sg.SerializableHdf5File(model_file, "w")
    with closing(output):
        if svm.save_serializable(output):
            print("Model saved: {}".format(model_file))
        else:
            print("Saving the model failed")
