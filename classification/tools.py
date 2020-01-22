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
            print("Model failed to load")
            exit(1)

    return svm


def performance_metrics(predict, labels, imbalance_rat):
    acc_measure = sg.AccuracyMeasure()
    acc = acc_measure.evaluate(predict, labels)
    prec = sg.PrecisionMeasure().evaluate(predict, labels)
    recall = sg.RecallMeasure().evaluate(predict, labels)

    metrics = [
        " -TP: {}".format(acc_measure.get_TP()),
        " -FP: {}".format(acc_measure.get_FP()),
        " -TN: {}".format(acc_measure.get_TN()),
        " -FN: {}".format(acc_measure.get_FN()),

        "Accuracy: {}".format(acc),
        "Precision: {}".format(prec),
        "Recall: {}".format(recall),
    ]

    def adjust_precision(acc_meas, imbalance_ratio):
        TP = acc_meas.get_TP()
        FP = acc_meas.get_FP()

        # Infer ratio of +/- classes in the validation set
        v = (TP + acc_meas.get_FN()) / (acc_meas.get_TN() + FP)
        # Ratio of +/- classes in real data
        r = imbalance_ratio

        return (r / v * TP) / (r / v * TP + (1 - r / 1 - v) * FP)

    if imbalance_rat:
        adj_prec = adjust_precision(acc_measure, imbalance_rat)
        metrics += [
            "Adjusted precision: {}".format(adj_prec),
            "Adj. precision calculated with imbalance ratio: {}".format(imbalance_rat)
        ]

    return metrics
