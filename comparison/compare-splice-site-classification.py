import pandas as pd


def main():
    """
    Compares splice site classification output of two models (usually NN and SVM, but can be 2 SVMs). Examines, how much
    classification is shared and whether the models tend to do similar mistakes.
    @nn_output  CSV output of neural network splice site classification. Columns [sequence, label, prediction]
    @svm_output CSV output of SVM splice site classification. Columns [sequence, label, prediction]
    @nn_adjustment  Switch to adjust NN windows to have the same length as SVM windows. Needed for joining of CSVs

    The files can be obtained via e.g.:
    python classify-splice-sites.py $FILE ../pipeline/bestmodels/ascomycota/model.hd5 70 70 acceptor \
                                -r 0.003  -c 12  > 'svm_acceptor_output.csv'
    :return: None
    """
    nn_output = f'nn_acceptor_output.csv'
    svm_output = f'svm_acceptor_output.csv'
    nn_adjustment = True

    nn_classification = pd.read_csv(nn_output, delimiter=';')
    svm_classification = pd.read_csv(svm_output, delimiter=';')

    if nn_adjustment:
        # Adjust range so the sequences of two files are equally long
        nn_classification['sequence'] = nn_classification[['sequence']].applymap(lambda seq: seq[128:270])

    nn_classification.drop_duplicates(inplace=True)
    svm_classification.drop_duplicates(inplace=True)

    joined_results = pd.merge(nn_classification, svm_classification, on='sequence')
    positives_df = joined_results[joined_results['label_x'] == 1]

    FNs_x = positives_df[positives_df['pred_x'] == -1].shape[0]
    FNs_y = positives_df[positives_df['pred_y'] == -1].shape[0]

    TPs_x = positives_df[positives_df['pred_x'] == 1].shape[0]
    TPs_y = positives_df[positives_df['pred_y'] == 1].shape[0]

    tp_tp_intersect = positives_df.query('pred_x == 1 & pred_y == 1').shape[0]
    print(f'True positives intersection {tp_tp_intersect}. '
          f'I.e {100 * tp_tp_intersect / TPs_x:.2f}% of model X '
          f'and {100 * tp_tp_intersect / TPs_y:.2f}% of model Y results')

    fn_fn_intersect = positives_df.query('pred_x == -1 & pred_y == -1').shape[0]
    print(f'False negatives intersection {fn_fn_intersect}. '
          f'The intersection makes {100 * fn_fn_intersect / FNs_x:.2f}% of model X false negatives '
          f'and {100 * fn_fn_intersect / FNs_y:.2f}% of model Y false negatives')

    fn_tp_intersect = positives_df.query('pred_x == -1 & pred_y == 1').shape[0]
    print(f'FN vs TP intersection {fn_tp_intersect}. '
          f'I.e {100 * fn_tp_intersect / FNs_x:.2f}% of Xs model FN is found by model Y ')

    tp_fn_intersect = positives_df.query('pred_x == 1 & pred_y == -1').shape[0]
    print(f'FN vs TP intersection {tp_fn_intersect}. '
          f'I.e {100 * tp_fn_intersect / FNs_y:.2f}% of Ys model FN is found by model X ')

    negatives_df = joined_results[joined_results['label_x'] == -1]
    FPs_x = negatives_df[negatives_df['pred_x'] == 1].shape[0]
    FPs_y = negatives_df[negatives_df['pred_y'] == 1].shape[0]

    fp_fp_intersect = negatives_df.query('pred_x == 1 & pred_y == 1').shape[0]
    print(f'False positives intersection {fp_fp_intersect}. '
          f'The intersection makes {100 * fp_fp_intersect / FPs_x:.2f}% of model X false positives '
          f'and {100 * fp_fp_intersect / FPs_y:.2f}% of model Y false positives')


if __name__ == '__main__':
    main()
