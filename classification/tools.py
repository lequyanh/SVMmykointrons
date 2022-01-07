from contextlib import closing

import pandas as pd
from sklearn.metrics import confusion_matrix


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
    import shogun as sg
    svm = sg.LibSVM()

    model_file = sg.SerializableHdf5File(filename, "r")

    with closing(model_file):
        if not svm.load_serializable(model_file):
            print("Model failed to load")
            exit(1)

    return svm


def split_data(data, test_size):
    #s = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=42)
    #for train_i, test_i in s.split(data.sequence, data.label):
    #    train, test = data.loc[train_i], data.loc[test_i]
    #    return train, test

    no_test_examples = int(test_size * data.shape[0])
    test_set = data[:no_test_examples]
    train_set = data[no_test_examples:]
    print(f'Positive samples test {sum(test_set["label"]==1)}')
    return train_set, test_set


def performance_metrics(labels, predictions, imbalance_rat):
    conf_matrix = confusion_matrix(labels, predictions)
    TN = conf_matrix[0][0]
    FN = conf_matrix[1][0]
    TP = conf_matrix[1][1]
    FP = conf_matrix[0][1]

    acc = (TP + TN) / (TP + FP + FN + TN)
    prec = TP / (TP + FP)
    recall = TP / (TP + FN)

    metrics = [
        " -TP: {}".format(TP),
        " -FP: {}".format(FP),
        " -TN: {}".format(TN),
        " -FN: {}".format(FN),

        "Accuracy: {}".format(acc),
        "Precision: {}".format(prec),
        "Recall: {}".format(recall)
    ]

    def adjust_precision(imbalance_ratio):
        # Infer ratio of +/- classes in the validation set
        v = (TP + FN) / (TN + FP)
        # Ratio of +/- classes in real data
        r = imbalance_ratio

        return (r / v * TP) / (r / v * TP + (1 - r / 1 - v) * FP)

    if imbalance_rat:
        adj_prec = adjust_precision(imbalance_rat)
        metrics += [
            "Adjusted precision: {}".format(adj_prec),
            "Adj. precision calculated with imbalance ratio: {}".format(imbalance_rat)
        ]

    return metrics
