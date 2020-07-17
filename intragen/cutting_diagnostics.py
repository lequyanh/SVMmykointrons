import logging

import pandas as pd
from pandas import DataFrame
from pandas.core.groupby import GroupBy

import config
import intragen_tools
from label_introns import get_introns_from_strand
from tools import performance_metrics

DONOR_SITE = 'donor'
ACCEPTOR_SITE = 'acceptor'


def main():
    fungi_name = 'Kocim1'
    model = 'SVM-intragen'
    strand = '+'
    folder = 'full_results'

    run_diagnostics(fungi_name, model, strand, folder)


def run_diagnostics(fungi_name: str, model: str, strand: str, folder: str = '.'):
    logging.basicConfig(
        level=logging.INFO,
        filename=f'{fungi_name}_intron_intragen_exploration.log',
        filemode='w'
    )

    exon_file = config.get_fungi_exons_positions(fungi_name)
    intron_file = config.get_fungi_intron_fasta(fungi_name)

    model = f'-{model}' if model else ''
    # Result or running the pipeline (intron-result.csv file)
    intron_annotations_file = f'{folder}/{fungi_name}_results/intron-result{model}.csv'
    # Result of running the pipeline (intron pruning step)
    cut_coords_file = f'{folder}/{fungi_name}_results/cut-coords{model}.csv'

    # NOTE - joined_df will be missing introns, which did not survive pairing phase
    joined_df = intragen_tools.join_on_position(cut_coords_file, intron_annotations_file)

    # Get the accuracy and recall of intron detection after the pruning step
    post_cut_accuracy_metrics(joined_df)

    # Determine intra-genic intron FP rate. We have to pass strand here as we don't know, against which exons to compare
    recall, exon_breaking_fpr = false_introns_exploration(joined_df, exon_file, intron_file, strand)

    return recall, exon_breaking_fpr


def post_cut_accuracy_metrics(joined: DataFrame) -> None:
    """
    Compares accuracy metrics of intron detection before and after cutting.
    The metrics may be different since there is a large number of overlaps and we need to choose which candidate to cut

    :param joined: DataFrame joined table of intron coordinates, their labels, their prediction and cut flag
    :return: None (only prints the metrics)
    """
    metrics_before = performance_metrics(joined['label'], joined['prediction'], None)
    metrics_before = '\n'.join(metrics_before)

    print('Before cut metrics')
    print(metrics_before)
    print('\n')

    metrics = performance_metrics(joined['label'], joined['cut'], None)
    metrics = '\n'.join(metrics)

    print('After cut metrics')
    print(metrics)
    print('\n')


def false_introns_exploration(
        joined: DataFrame,
        exon_file: str,
        intron_file: str,
        strand: str
) -> (float, float):
    """
    Explores, where false introns fall into. Determines what portion of them lie inside exons.
    Prints adjusted FP intron rate, where only FP inside exons are considered.

    :param intron_file: Filename of fungi intron sequences
    :param joined: DataFrame joined table of intron coordinates, their labels, their prediction and cut flag
    :param exon_file: Filename of fungi exon positions
    :param strand: +/- so we know, which exons to pick
    :return: None (only prints the metrics)
    """
    exon_pos_df = pd.read_csv(exon_file, delimiter=';')
    if 'strand' in exon_pos_df.columns:
        # For legacy reasons - some exon files don't yet have a "strand" column
        exon_pos_df = exon_pos_df[exon_pos_df['strand'] == strand]

    exon_scaff_grouped = exon_pos_df.groupby(by='scaffold')

    # See, where potential mistakes of pruning can be (intron dataset before classification)
    potential_fp_candidates = joined.query('label == -1')
    potential_fp_grouped = potential_fp_candidates.groupby(by='scaffold')
    print(f'--------------------------- DATASET ----------------------------')
    no_cuts = intraexonic_cuts_count(exon_scaff_grouped, potential_fp_grouped)
    print(
        f'{no_cuts}/{len(potential_fp_candidates)} false CANDIDATES are inside exon => '
        f'the ratio of exon-breaking candidates in dataset is {100 * no_cuts / len(potential_fp_candidates):.2f}%.\n'
    )

    # See, where potential mistakes of pruning can be (intron classification intra-genic FP)
    labeled_true_df = joined.query('prediction == 1')
    intron_classification_fps = labeled_true_df.query('label == -1')
    fps_grouped = intron_classification_fps.groupby(by='scaffold')
    print(f'-------------------------- CLASSIFICATION --------------------------')
    no_cuts = intraexonic_cuts_count(exon_scaff_grouped, fps_grouped)
    print(
        f'{no_cuts}/{len(intron_classification_fps)} false POSITIVES are inside exon => '
        f'the proportion of exon-breaking FP among all FP is {no_cuts / len(intron_classification_fps):.2f}.\n'
        f'Exon-breaking FPR is {100 * no_cuts / labeled_true_df.shape[0]:.2f}%'
    )

    # See, where false intron cuts happened
    all_cuts_count = joined.query('cut == 1').shape[0]
    false_positive_cuts_df = intron_classification_fps.query('cut == 1')
    fp_cuts_grouped = false_positive_cuts_df.groupby(by='scaffold')

    print(f'-------------------------- CUTTING -------------------------------------')
    no_cuts = intraexonic_cuts_count(exon_scaff_grouped, fp_cuts_grouped)
    print(
        f'{no_cuts}/{len(false_positive_cuts_df)} false CUTS => '
        f'the proportion of exon-breaking cuts among all false cuts is {no_cuts / len(false_positive_cuts_df):.2f}.\n'
        f'After cut exon-breaking FPR is {100 * no_cuts / all_cuts_count:.2f}% as there are {all_cuts_count} total cuts'
    )

    print(f'--------------------------------Recall and exon-breaking precision ---------------------------------------')
    true_cuts = joined.query('cut == 1 and label == 1').shape[0]
    true_all = get_introns_from_strand(intron_file, strand='+')

    print(f'Correctly cut {true_cuts}/{len(true_all)} ({strand} strand) introns.\n'
          f'Interfered with {no_cuts}/{exon_pos_df.shape[0]} ({strand} strand) exons.\n'
          f'Ratio {true_cuts / no_cuts:.2f}')

    recall = true_cuts / len(true_all) if true_cuts > 0 else 0
    exon_breaking_fpr = no_cuts / exon_pos_df.shape[0]

    return recall, exon_breaking_fpr


def intraexonic_cuts_count(exon_grouped: GroupBy, cuts_grouped: GroupBy) -> int:
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
            logging.warning(f'No exons in scaffold {scaffold}, or scaffold missing')
            continue

        exon_cuts, exon_length = intragen_tools.intraexonic_cuts_in_scaffold(scaff_exon_positions, cut_positions)

        exon_cuts_no = sum(exon_cuts)
        logging.info(f'{scaffold}: {exon_cuts_no} within exon cuts. '
                     f'This makes {exon_cuts_no * 1000 / exon_length:.2f} cuts per 1k exon bp')

        total_intraexon_cuts += exon_cuts_no

    return total_intraexon_cuts


if __name__ == "__main__":
    main()
