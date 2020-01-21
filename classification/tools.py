import pandas as pd
import shogun as sg
from contextlib import closing


def read_data(filename, window):
    # center - lwin = start of the window (inclusively)
    # center + rwin + 2 = end of the window (exclusively)
    lwin = window[0]
    rwin = window[1]

    data = pd.read_csv(filename, sep=';')
    center = int(len(data.iloc[0]['sequence']) / 2 - 1)  # e.g. 402 / 2 - 1 = 200
    data.sequence = data.sequence.str[center - lwin: center + rwin + 2]

    return data


def read_model(filename):
    svm = sg.LibSVM()

    model_file = sg.SerializableHdf5File(filename, "r")

    with closing(model_file):
        if not svm.load_serializable(model_file):
            exit(1)

    return svm
