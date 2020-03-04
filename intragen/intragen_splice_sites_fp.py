import logging

import intragen_tools

NEWSEQUENCES_LOC = '/home/anhvu/Desktop/mykointrons-data/new-sequences'
DONOR_SITE = 'donor'
ACCEPTOR_SITE = 'acceptor'


def main():
    shroom_name = 'Kocim1'
    model = 'SVM'

    logging.basicConfig(
        level=logging.INFO,
        filename=f'{shroom_name}_splice-site_intragen_exploration.log',
        filemode='w'
    )

    exon_csv = f'{NEWSEQUENCES_LOC}/{shroom_name}/{shroom_name}_exon_positions.csv'
    introns_fasta = f'{NEWSEQUENCES_LOC}/{shroom_name}/{shroom_name}-introns.fasta'

    donor_dataset_csv = f'{shroom_name}/{shroom_name}-splice-site-donor-dataset.csv'
    donor_result_csv = f'{shroom_name}/{shroom_name}-splice-site-donor-result-{model}.csv'

    donor_intragen_fps = intragenic_splicesite_fpr(
        donor_dataset_csv,
        donor_result_csv,
        introns_fasta,
        exon_csv,
        DONOR_SITE
    )

    acc_dataset_csv = f'{shroom_name}/{shroom_name}-splice-site-acceptor-dataset.csv'
    acc_result_csv = f'{shroom_name}/{shroom_name}-splice-site-acceptor-result-{model}.csv'

    acceptor_intragen_fps = intragenic_splicesite_fpr(
        acc_dataset_csv,
        acc_result_csv,
        introns_fasta,
        exon_csv,
        ACCEPTOR_SITE
    )


def intragenic_splicesite_fpr(
        classification_data_csv,
        classification_result_csv,
        introns_fasta,
        exon_csv,
        site: str
) -> float:
    """
    Determine intra-genic FP ratio for given splice site.
    Does so by merging classification data (windows and their positions) with their classification and labels.
    Lastly, extends the merged DataFrame with information, whether the candidate lies inside exons or not.
    """
    merged_df = intragen_tools.merge_classification_data(
        classification_data_csv,
        classification_result_csv,
        introns_fasta,
        exon_csv,
        site
    )

    predict_trues = merged_df.query('prediction == 1')
    total_fps = predict_trues.query('label == -1')
    intra_exon_fps = total_fps.query('in_exon == 1').shape[0]

    fpr = 100 * total_fps.shape[0] / predict_trues.shape[0]
    intragen_fpr = 100 * intra_exon_fps / predict_trues.shape[0]

    print(f'{site} site predicted trues {predict_trues.shape[0]},\n'
          f'from which false positives {total_fps.shape[0]} (FP rate {fpr:.2f}%),\n'
          f'from which intra-genic FPs {intra_exon_fps} (intra-genic FP ratio {intragen_fpr:.2f}%)\n')

    return intragen_fpr


if __name__ == '__main__':
    main()
