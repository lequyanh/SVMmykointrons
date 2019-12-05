import numpy as np
import pandas as pd
import sys
from contextlib import closing
import shogun as sg


def read_data(filename, window):
    # center - lwin = start of the window (inclusively)
    # center + rwin + 2 = end of the window (exclusively)
    center = 150
    lwin = window[0]
    rwin = window[1]

    data = pd.read_csv(filename, sep=';')
    data.sequence = data.sequence.str[center - lwin: center + rwin + 2]

    #    print("Loaded {} samples".format(len(data)))
    #    print("Using sequence window [{}, {}]".format(lwin, rwin))

    return data


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
    print("Accuracy: {}".format(acc.evaluate(predict, labels)))
    print(" -TP: {}".format(acc.get_TP()))
    print(" -FP: {}".format(acc.get_FP()))
    print(" -TN: {}".format(acc.get_TN()))
    print(" -FN: {}".format(acc.get_FN()))

    output = sg.SerializableHdf5File(model_file, "w")
    with closing(output):
        if svm.save_serializable(output):
            print("Model saved: {}".format(model_file))
        else:
            print("Saving the model failed")
