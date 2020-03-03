import logging

import pandas as pd
from pandas import DataFrame

import intragen_tools
from extract_tools import true_donorac_positions

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
    merged_df = merge_classification_data(
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


def merge_classification_data(
        classification_data_csv,
        classification_result_csv,
        introns_fasta,
        exon_csv,
        site: str
) -> DataFrame:
    classification_data_df = pd.read_csv(classification_data_csv, delimiter=';')
    # Remove unnecessary sequence column
    classification_data_df.drop(columns='sequence', inplace=True)

    merged_df = merge_data_prediction_labels(
        classification_data_df,
        classification_result_csv,
        introns_fasta,
        site
    )

    merged_df.rename(columns={"position": "start"}, inplace=True)  # Rename so it matches required interface
    merged_df = add_in_exon_column(merged_df, exon_csv)

    return merged_df


def merge_data_prediction_labels(
        classification_data_df: DataFrame,
        classification_result_csv: str,
        introns_fasta: str,
        site: str
) -> DataFrame:
    """
    Merges donor/acceptor dataset with predictions and labels.
    The goal is to have attributes 'scaffold', 'position', 'prediction', 'label' for each candidate window

    :param classification_data_df: DataFrame of candidate splice site (windows) ['scaffold', 'position']
    :param classification_result_csv: CSV file with positively classified splice site positions ['scaffold', 'position']
    :param introns_fasta: FASTA file with introns
    :param site: Donor/Acceptor
    :return: Merged DataFrame with ['scaffold', 'position', 'prediction', 'label']
    """
    classification_result_df = pd.read_csv(classification_result_csv, delimiter=';')
    # Add new column and set it to be ones
    classification_result_df['prediction'] = 1

    # Create a DF with true donor/acceptor positions (i.e. DataFrame of candidate labels)
    true_donor_pos, true_acceptor_pos = true_donorac_positions(introns_fasta)

    if site == DONOR_SITE:
        true_splicesite_df = intragen_tools.splice_site_positions_dict_to_df(true_donor_pos)
    else:
        true_splicesite_df = intragen_tools.splice_site_positions_dict_to_df(true_acceptor_pos)
    # Add new column and set it to be ones
    true_splicesite_df['label'] = 1

    # Merge all the 3 data sets (input dataset, classification and labels) together
    merged_df = classification_data_df.merge(true_splicesite_df, how="left").fillna(-1)
    merged_df = merged_df.merge(classification_result_df, how="left").fillna(-1)

    return merged_df


def add_in_exon_column(
        merged_df: DataFrame,
        exon_csv: str
) -> DataFrame:
    """
    Extends a given DataFrame of candidate splice sites with column 'in_exon'.

    :param merged_df: DataFrame with columns ['scaffold', 'start'], preferably also with ['prediction', 'label']
    :param exon_csv: CSV with exon positions
    :return: Merged DataFrame with ['scaffold', 'start', 'prediction', 'label', 'in_exon']
    """
    merged_df['in_exon'] = -1  # New column of "in_exon" indicators. Initialize to -1 and update those inside exons

    # Create DataFrame with exon positions. The DF will be used to determine intra-exonic false predictions
    exon_pos_df = pd.read_csv(exon_csv, delimiter=';')

    exon_grouped = exon_pos_df.groupby(by='scaffold')
    merged_grouped = merged_df.groupby(by='scaffold')

    for scaffold, positions in merged_grouped:
        try:
            scaff_exon_positions = exon_grouped.get_group(scaffold)
        except KeyError:
            logging.warning(f'No exons in scaffold {scaffold}, or scaffold missing')
            continue

        # Pick positions (within given scaffold) of false positives and find, whether they fall into exon or not
        filter_false_positives = (merged_df.scaffold == scaffold) & (merged_df.label == -1) & (
                merged_df.prediction == 1)
        cut_positions = merged_df[filter_false_positives]
        exon_cut_indicators, _ = intragen_tools.intraexonic_cuts_in_scaffold(scaff_exon_positions, cut_positions)

        merged_df.loc[filter_false_positives, 'in_exon'] = list(
            map(lambda x: -1 if x is False else 1, exon_cut_indicators)
        )

    return merged_df


if __name__ == '__main__':
    main()
