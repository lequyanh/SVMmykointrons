import itertools
from bisect import bisect
from typing import Tuple, List

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.core.groupby import GroupBy

from extract_tools import true_donorac_positions
from tools import performance_metrics

NEWSEQUENCES_LOC = '/home/anhvu/Desktop/mykointrons-data/new-sequences'


def main():
    """
    @cut_coords_file  Result of intron pruning. While pruning, record coords which are cut. Can be obtained by running
                            echo "start;end" > cut_coords.csv &
                            python ../pruning/prune_probabilistic.py
                                    ~/Desktop/mykointrons-data/data/Assembly/Kocim1_AssemblyScaffolds.fasta
                                    intron-result.csv
                                    ../statistics/intron-lenghts-statistics/basidiomycota-intron-lens.txt
                            >> cut_coords.csv
    @intron_annotation_file     Result or running pipeline with an exception - do not discard negatively labeled
                                introns candidates
    """
    shroom_name = 'Thega1'
    exon_file = f'{NEWSEQUENCES_LOC}/{shroom_name}/{shroom_name}_exon_positions.csv'

    intron_annotations_file = f'{shroom_name}_intron-result-NN.csv'
    cut_coords_file = f'{shroom_name}_cut_coords-NN.csv'

    # Get the accuracy and recall of intron detection after the pruning step
    joined_df = join_on_position(cut_coords_file, intron_annotations_file)
    get_post_cut_accuracy_metrics(joined_df)

    exon_pos_df = pd.read_csv(exon_file, delimiter=';')
    exon_scaff_grouped = exon_pos_df.groupby(by='scaffold')

    # See, where potential mistakes of pruning can be (intron dataset before classification)
    negative_intron_candidates = joined_df[joined_df.label == -1]
    intron_candidates_grouped = negative_intron_candidates.groupby(by='scaffold')
    print(f'--------------------- All false potential intra-exon cuts ---------------------------')
    no_cuts = determine_num_intraexonic_cuts(exon_scaff_grouped, intron_candidates_grouped)
    print(
        f'>> Cuts: {no_cuts} out of {len(negative_intron_candidates)} ({no_cuts / len(negative_intron_candidates):.2f})\n')

    # See, where false intron cuts happened
    false_positive_cuts_df = joined_df.query('pred == 1 and label == -1 and cut == 1')
    fp_cuts_grouped = false_positive_cuts_df.groupby(by='scaffold')

    print(f'--------------------- False positive intra-exon cuts ---------------------------')
    no_cuts = determine_num_intraexonic_cuts(exon_scaff_grouped, fp_cuts_grouped)
    print(f'>> Cuts: {no_cuts} out of {len(false_positive_cuts_df)} ({no_cuts / len(false_positive_cuts_df):.2f})')


def main2():
    shroom_name = "Kocim1"
    exon_file = f'{NEWSEQUENCES_LOC}/{shroom_name}/{shroom_name}_exon_positions.csv'
    introns_fasta = f'{NEWSEQUENCES_LOC}/{shroom_name}/{shroom_name}-introns.fasta'

    classification_dataset_f = f'{shroom_name}-splice-site-donor-dataset.csv'
    classification_result_f = f'{shroom_name}-splice-site-donor-result-svm.csv'

    classification_data_df = pd.read_csv(classification_dataset_f, delimiter=';')
    # Remove unnecessary sequence column and merge the classification dataset with labels and predictions
    classification_data_df.drop(columns='sequence', inplace=True)

    classification_result_df = pd.read_csv(classification_result_f, delimiter=';')
    # Add new column and set it to be ones
    classification_result_df['prediction'] = 1

    # Create a DF with true donor/acceptor positions (i.e. DataFrame of candidate labels)
    true_donor_pos, true_acceptor_pos = true_donorac_positions(introns_fasta)

    true_donors_df = splice_site_positions_dict_to_df(true_donor_pos)
    true_donors_df['label'] = 1  # Add new column and set it to be ones

    # Merge all the 3 data sets (input dataset, classification and labels) together
    merged_df = classification_data_df.merge(true_donors_df, how="left").fillna(-1)
    merged_df = merged_df.merge(classification_result_df, how="left").fillna(-1)

    merged_df.rename(columns={"position": "start"}, inplace=True)  # Rename so it matches required interface
    merged_df['in_exon'] = -1  # New column of "in_exon" indicators. Initialize to -1 and update those inside exons

    # Create DataFrame with exon positions. The DF will be used to determine intra-exonic false predictions
    exon_pos_df = pd.read_csv(exon_file, delimiter=';')

    exon_grouped = exon_pos_df.groupby(by='scaffold')
    merged_grouped = merged_df.groupby(by='scaffold')

    for scaffold, positions in merged_grouped:
        try:
            scaff_exon_positions = exon_grouped.get_group(scaffold)
        except KeyError:
            print(f'No exons in scaffold {scaffold}, or scaffold missing')
            continue

        # Pick positions (within given scaffold) of false positives and find, whether they fall into exon or not
        filter_false_positives = (merged_df.scaffold == scaffold) & (merged_df.label == -1) & (
                    merged_df.prediction == 1)
        cut_positions = merged_df[filter_false_positives]
        exon_cut_indicators, _ = get_intraexonic_cuts(scaff_exon_positions, cut_positions)

        merged_df.loc[filter_false_positives, 'in_exon'] = list(
            map(lambda x: -1 if x is False else 1, exon_cut_indicators)
        )


def splice_site_positions_dict_to_df(splicesite_positions: dict):
    per_scaffold_dfs = []
    for scaffold, positions in splicesite_positions.items():
        data = {'scaffold': [scaffold] * len(positions), 'position': positions}
        per_scaffold_dfs.append(pd.DataFrame(data))

    return pd.concat(per_scaffold_dfs)


def determine_num_intraexonic_cuts(exon_grouped: GroupBy, cuts_grouped: GroupBy) -> int:
    """
    Determine, how many (potential) cuts happened inside exon sequences.
    :param exon_grouped: GroupBy exon position for each scaffold
    :param cuts_grouped: GroupBy cut positions for each scaffold
    :return: Number of (potential) cuts inside exons
    """
    total_intraexon_cuts = 0
    for scaffold, cut_positions in cuts_grouped:
        try:
            scaff_exon_positions = exon_grouped.get_group(scaffold)
        except KeyError:
            print(f'No exons in scaffold {scaffold}, or scaffold missing')
            continue

        exon_cuts, exon_length = get_intraexonic_cuts(scaff_exon_positions, cut_positions)

        exon_cuts_no = sum(exon_cuts)
        print(f'{scaffold}: {exon_cuts_no} within exon cuts. '
              f'This makes {exon_cuts_no * 1000 / exon_length:.2f} cuts per 1k exon bp')

        total_intraexon_cuts += exon_cuts_no

    return total_intraexon_cuts


def get_intraexonic_cuts(scaff_exon_positions: DataFrame, query_positions: DataFrame) -> Tuple[List[int], int]:
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


def get_post_cut_accuracy_metrics(joined: DataFrame) -> None:
    """
    Compares accuracy metrics of intron detection before and after cutting.
    The metrics may be different since there is a large number of overlaps and we need to choose which candidate to cut
    :param joined: DataFrame joined table of intron coordinates, their labels, their prediction and cut flag
    :return: None (only prints the metrics)
    """
    metrics_before = performance_metrics(joined['label'], joined['pred'], None)
    metrics_before = '\n'.join(metrics_before)

    print('Before cut metrics')
    print(metrics_before)
    print('\n')

    metrics = performance_metrics(joined['label'], joined['cut'], None)
    metrics = '\n'.join(metrics)

    print('After cut metrics')
    print(metrics)


def join_on_position(cut_coords_file: str, intron_annotations_file: str) -> DataFrame:
    """
    Joins a table intron cuts (their coordinates) and a table of labeled intron candidates.
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


if __name__ == "__main__":
    main2()
