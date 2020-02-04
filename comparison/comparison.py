import pandas as pd

nn_output = f'nn_acceptor_output.csv'
svm_output = f'svm_acceptor_output.csv'

nn_classification = pd.read_csv(nn_output, delimiter=';')
svm_classification = pd.read_csv(svm_output, delimiter=';')

# Adjust range so the sequences of two files are equally long
nn_classification['sequence'] = nn_classification[['sequence']].applymap(lambda seq: seq[128:270])

intersect_whole = pd.merge(nn_classification, svm_classification, how='inner')
print(f'Intersection in classification: {len(intersect_whole)} out of {len(nn_classification)}.'
      f'Ratio {100 * len(intersect_whole) / len(nn_classification)}%')


def get_intersection(filter_nn, filter_svm):
    nn_view = nn_classification[filter_nn]
    svm_view = svm_classification[filter_svm]

    # Drop predictions so the join makes sense
    nn_view = nn_view.drop(columns='pred')
    svm_view = svm_view.drop(columns='pred')

    intersect = pd.merge(nn_view, svm_view, how='inner')

    return intersect, nn_view, svm_view


filter_nn_FN = (nn_classification['pred'] == -1) & (nn_classification['label'] == 1)
filter_svm_TP = (svm_classification['pred'] == 1) & (nn_classification['label'] == 1)

intersect_nnFN_svmTP, nn_fn, _ = get_intersection(filter_nn_FN, filter_svm_TP)
print(f'Intersection between SVM TP and NN FN is {len(intersect_nnFN_svmTP)}/{len(nn_fn)}. '
      f'I.e {100 * len(intersect_nnFN_svmTP) / len(nn_fn)}% of NN false negatives are detected by SVM')

# -------------------------------------------------------------------------------------------------------
filter_nn_TP = (nn_classification['pred'] == 1) & (nn_classification['label'] == 1)
filter_svm_FN = (svm_classification['pred'] == -1) & (nn_classification['label'] == 1)

intersect_svmFN_nnTP, _, svm_fn = get_intersection(filter_nn_TP, filter_svm_FN)
print(f'Intersection between SVM TP and NN FN is {len(intersect_svmFN_nnTP)}/{len(svm_fn)}. '
      f'I.e {100 * len(intersect_svmFN_nnTP) / len(svm_fn)}% of SVM false negatives are detected by MM')

# ------------------------------------------SHARING IS CARING-----------------------------------------------------------
filter_nn_FP = (nn_classification['pred'] == 1) & (nn_classification['label'] == -1)
filter_svm_FP = (svm_classification['pred'] == 1) & (nn_classification['label'] == -1)

intersect_svmFP_nnFP, nn_fp, svm_fp = get_intersection(filter_nn_FP, filter_svm_FP)
print(f'Intersection between SVM FP and NN FP is {len(intersect_svmFP_nnFP)}. '
      f'I.e {100 * len(intersect_svmFP_nnFP) / len(svm_fp)}% of SVM false positives are shared with NN '
      f'and {100 * len(intersect_svmFP_nnFP) / len(nn_fp)}% of NN false positives are shared with SVM')

intersect_svmTP_nnTP, nn_tp, svm_tp = get_intersection(filter_nn_TP, filter_svm_TP)
print(f'Intersection between SVM TP and NN TP is {len(intersect_svmTP_nnTP)}. '
      f'I.e {100 * len(intersect_svmTP_nnTP) / len(svm_tp)}% of SVM true positives are shared with NN '
      f'and {100 * len(intersect_svmTP_nnTP) / len(nn_tp)}% of NN true positives are shared with SVM')

intersect_svmFN_nnFN, nn_fn, svm_fn = get_intersection(filter_nn_FN, filter_svm_FN)
print(f'Intersection between SVM FN and NN FN is {len(intersect_svmFN_nnFN)}. '
      f'I.e {100 * len(intersect_svmFN_nnFN) / len(svm_fn)}% of SVM false negatives are shared with NN '
      f'and {100 * len(intersect_svmFN_nnFN) / len(nn_fn)}% of NN false negatives are shared with SVM')

