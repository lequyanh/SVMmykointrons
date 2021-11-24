import itertools
import logging
from bisect import bisect
from typing import Tuple, List

import numpy as np
import pandas as pd
from pandas import DataFrame

from extract_tools import true_donorac_positions

DONOR_SITE = 'donor'
ACCEPTOR_SITE = 'acceptor'


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

    merged_df = _merge_data_prediction_labels(
        classification_data_df,
        classification_result_csv,
        introns_fasta,
        site
    )

    merged_df.rename(columns={"position": "start"}, inplace=True)  # Rename so it matches required interface
    merged_df = _add_in_exon_column(merged_df, exon_csv)

    return merged_df


def _add_in_exon_column(
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
        exon_cut_indicators, _ = intraexonic_cuts_in_scaffold(scaff_exon_positions, cut_positions)

        merged_df.loc[filter_false_positives, 'in_exon'] = list(
            map(lambda x: -1 if x is False else 1, exon_cut_indicators)
        )

    return merged_df


def _merge_data_prediction_labels(
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
        true_splicesite_df = splice_site_positions_dict_to_df(true_donor_pos)
    else:
        true_splicesite_df = splice_site_positions_dict_to_df(true_acceptor_pos)
    # Add new column and set it to be ones
    true_splicesite_df['label'] = 1

    # Merge all the 3 data sets (input dataset, classification and labels) together
    merged_df = classification_data_df.merge(true_splicesite_df, how="left").fillna(-1)
    merged_df = merged_df.merge(classification_result_df, how="left").fillna(-1)

    return merged_df


def intraexonic_cuts_in_scaffold(scaff_exon_positions: DataFrame, query_positions: DataFrame) -> Tuple[List[int], int]:
    """
    Determines which of @query_positions lie inside an exon (within a particular scaffold)
    :param scaff_exon_positions: Exon positions in a single scaffold (SINGLE)
    :param query_positions: Positions we want to evaluate
    :return: List of True/False indicators, the length of exon
    """
    start_end_tuplelist = zip(scaff_exon_positions.start, scaff_exon_positions.end)
    start_end_tuplelist = sorted(start_end_tuplelist, key=lambda tup: tup[0])

    start_end_flattened = list(itertools.chain(*start_end_tuplelist))
    exon_length = sum([se[1] - se[0] for se in start_end_tuplelist])

    def inside_exon_cut(cut_start) -> int:
        i = bisect(start_end_flattened, cut_start)
        return i % 2 == 1

    cut_starts = query_positions['start']
    exon_cuts = [inside_exon_cut(cut_start) for cut_start in cut_starts]

    return exon_cuts, exon_length


def join_on_position(cut_coords_file: str, intron_annotations_file: str) -> DataFrame:
    """
    Joins a table of intron cuts (their coordinates) and a table of labeled intron candidates.
    The goal is to be able to see for each candidate his label, prediction and cut decision.
    :param cut_coords_file: CSV file with coordinates of performed cuts
    :param intron_annotations_file: CSV file of labeled intron candidates as well their prediction flags
    :return: DataFrame with columns [scaffold, start, end, label, prediction, cut]
    """
    cut_coords_df = pd.read_csv(cut_coords_file, delimiter=';')
    intron_coords_df = pd.read_csv(intron_annotations_file, delimiter=';')

    cut_coords_df['cut'] = np.ones(len(cut_coords_df))

    joined = intron_coords_df.merge(cut_coords_df, how="left")
    joined[['cut']] = joined[['cut']].fillna(value=-1)

    return joined


def splice_site_positions_dict_to_df(splicesite_positions: dict) -> DataFrame:
    """
    Convert dictionary in form of {scaffold: [123, 12, 34], ...} to a data frame.
    :param splicesite_positions: Dictionary to convert
    :return: DataFrame with 'scaffold' column
    """
    per_scaffold_dfs = []
    for scaffold, positions in splicesite_positions.items():
        data = {'scaffold': [scaffold] * len(positions), 'position': positions}
        per_scaffold_dfs.append(pd.DataFrame(data))

    return pd.concat(per_scaffold_dfs)