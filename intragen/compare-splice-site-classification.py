import pandas as pd

import config
import intragen_tools

DONOR_SITE = 'donor'
ACCEPTOR_SITE = 'acceptor'


def main():
    """
    Compares splice site classification output of two models (usually NN and SVM, but can be 2 SVMs). Examines, how much
    classification is shared and whether the models tend to do similar mistakes.
    :return:
    """
    fungi_name = 'Thega1'
    site = 'donor'

    model1 = 'SVM-intragen'
    model2 = 'NN'

    exon_csv = config.get_fungi_exons_positions(fungi_name)
    introns_fasta = config.get_fungi_intron_fasta(fungi_name)

    dataset_csv = f'../intragen/{fungi_name}/{fungi_name}-splice-site-{site}-dataset.csv'
    model1_result_csv = f'../intragen/{fungi_name}/{fungi_name}-splice-site-{site}-result-{model1}.csv'
    model2_result_csv = f'../intragen/{fungi_name}/{fungi_name}-splice-site-{site}-result-{model2}.csv'

    merged_df1 = intragen_tools.merge_classification_data(
        dataset_csv,
        model1_result_csv,
        introns_fasta,
        exon_csv,
        DONOR_SITE
    )

    merged_df2 = intragen_tools.merge_classification_data(
        dataset_csv,
        model2_result_csv,
        introns_fasta,
        exon_csv,
        DONOR_SITE
    )

    joined_results = pd.merge(merged_df1, merged_df2, on=['scaffold', 'start', 'label'])
    positives_df = joined_results[joined_results['label'] == 1]
    TPs_x = positives_df[positives_df['prediction_x'] == 1].shape[0]
    TPs_y = positives_df[positives_df['prediction_y'] == 1].shape[0]

    tp_tp_intersect = positives_df.query('prediction_x == 1 & prediction_y == 1').shape[0]
    print(f'True positives intersection {tp_tp_intersect}. '
          f'I.e {100 * tp_tp_intersect / TPs_x:.2f}% of {model1} model '
          f'and {100 * tp_tp_intersect / TPs_y:.2f}% of {model2} model results')

    negatives_df = joined_results[joined_results['label'] == -1]
    FPs_x = negatives_df[negatives_df['prediction_x'] == 1].shape[0]
    FPs_y = negatives_df[negatives_df['prediction_y'] == 1].shape[0]

    fp_fp_intersect = negatives_df.query('prediction_x == 1 & prediction_y == 1').shape[0]
    print(f'False positives intersection {fp_fp_intersect}. '
          f'The intersection makes {100 * fp_fp_intersect / FPs_x:.2f}% of {model1} model false positives '
          f'and {100 * fp_fp_intersect / FPs_y:.2f}% of {model2} model false positives')

    FPs_x = negatives_df.query('prediction_x == 1 and in_exon_x == 1').shape[0]
    FPs_y = negatives_df.query('prediction_y == 1 and in_exon_y == 1').shape[0]

    fp_fp_intersect = negatives_df.query(
        'prediction_x == 1 & prediction_y == 1 & in_exon_x == 1 & in_exon_y == 1'
    ).shape[0]
    print(f'Intragenic false positives intersection {fp_fp_intersect}. '
          f'The intersection makes {100 * fp_fp_intersect / FPs_x:.2f}% of {model1} model false positives '
          f'and {100 * fp_fp_intersect / FPs_y:.2f}% of {model2} model false positives')


if __name__ == '__main__':
    main()
