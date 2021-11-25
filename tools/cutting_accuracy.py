import logging
import sys

import pandas as pd
from Bio import SeqRecord, SeqIO
from pandas import DataFrame
from pandas.core.groupby import GroupBy

import intragen_tools


def main():
    cut_coords_file = sys.argv[1]
    intron_file = sys.argv[2]
    exon_file = sys.argv[3]
    strand = sys.argv[4]
    assembly = sys.argv[5] if len(sys.argv) == 6 else None

    if strand not in ['-', '+']:
        strand = '-' if strand == 'minus' else '+'

    run_diagnostics(cut_coords_file, intron_file, exon_file, strand, assembly)


def run_diagnostics(cut_coords_file: str, intron_file: str, exon_file, strand: str, assembly_fasta: str):
    logging.basicConfig(
        level=logging.INFO,
        filename=f'cut_stats.log',
        filemode='w'
    )

    # Result of running the pipeline (intron pruning step)
    cut_introns_df = pd.read_csv(cut_coords_file, delimiter=';')
    cut_introns_df['cut'] = 1

    # Convert cut coordinates if strand is minus
    if strand == '-':
        with open(assembly_fasta) as f:
            lens = {sr.id: len(sr.seq) for sr in SeqIO.parse(f, 'fasta')}
        # Create a DF from with scaffold - scaffold_length pairs
        #      scaffold   length
        #    scaffold_1  2013504
        scaf_lengths = pd.DataFrame \
            .from_dict(lens, orient='index', columns=['length']) \
            .reset_index() \
            .rename(columns={'index': 'scaffold'})
        # Merge cuts with scaffold lengths and convert reverse strand coordinates
        cut_length_df = pd.merge(cut_introns_df, scaf_lengths, on='scaffold')
        cut_introns_df['start'] = cut_length_df['length'] - cut_length_df['end'] + 1
        cut_introns_df['end'] = cut_length_df['length'] - cut_length_df['start'] + 1

    # Load true intron locations
    true_introns_df = get_introns_from_strand(intron_file, strand)

    joined_df = pd.merge(true_introns_df, cut_introns_df,
                         on=['scaffold', 'start', 'end'],
                         how="outer")
    joined_df[['label']] = joined_df[['label']].fillna(value=-1)
    joined_df[['cut']] = joined_df[['cut']].fillna(value=-1)

    # Get the accuracy and recall of intron detection after the pruning step
    # post_cut_accuracy_metrics(joined_df)

    # Determine intra-genic intron FP rate. We have to pass strand here as we don't know, against which exons to compare
    true_cuts, exon_cuts, all_cuts, all_introns, detectable_introns, exons = \
        false_introns_exploration(joined_df, exon_file, strand)

    return true_cuts, exon_cuts, all_cuts, all_introns, detectable_introns, exons


def get_introns_from_strand(true_introns_fasta: str, strand: str):
    def is_strand(seq_rec: SeqRecord, sign: str):
        return seq_rec.description.split(' ')[1] == sign

    with open(true_introns_fasta, 'r') as intron_fasta_f:
        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))
        print(f'Parsing introns from {true_introns_fasta} file.')
        print(f'Total number of introns: {len(introns_seqrecords)}')

        introns = list(filter(lambda i_rec: is_strand(i_rec, strand), introns_seqrecords))
        introns_df = pd.DataFrame({
            'scaffold': [intron.description.split(' ')[0] for intron in introns],
            'start': [int(intron.description.split(' ')[2]) for intron in introns],
            'end': [int(intron.description.split(' ')[3]) for intron in introns],
            'label': 1
        })

        print(f'Number of true introns on the given ({strand}) strand: {len(introns)}')

        return introns_df


# def post_cut_accuracy_metrics(joined: DataFrame) -> None:
#     """
#     Compares accuracy metrics of intron detection before and after cutting.
#     The metrics may be different since "there is a large number of overlaps and we need to choose which candidate to cut
#     :param joined: DataFrame joined table of intron coordinates, their labels, their prediction and cut flag
#     :return: None (only prints the metrics)
#     """
#     metrics_before = performance_metrics(joined['label'], joined['prediction'], None)
#     metrics_before = '\n'.join(metrics_before)
#
#     print('Before cut metrics')
#     print(metrics_before)
#     print('\n')
#
#     metrics = performance_metrics(joined['label'], joined['cut'], None)
#     metrics = '\n'.join(metrics)
#
#     print('After cut metrics')
#     print(metrics)
#     print('\n')


def false_introns_exploration(
        joined: DataFrame,
        exon_file: str,
        strand: str
) -> (float, float):
    """
    Explores, where false introns fall into. Determines what portion of them lie inside exons.
    Prints adjusted FP intron rate, where only FP inside exons are considered.
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

    no_cuts = intraexonic_cuts_count(exon_scaff_grouped, potential_fp_grouped)
    print(
        f'{no_cuts}/{len(potential_fp_candidates)} false CANDIDATES are inside exon => '
        f'the ratio of exon-breaking candidates in dataset is {100 * no_cuts / len(potential_fp_candidates):.2f}%.\n'
    )

    # See, where false intron cuts happened
    all_cuts_count = joined.query('cut == 1').shape[0]
    false_positive_cuts_df = joined.query('cut == 1 and label == -1')
    fp_cuts_grouped = false_positive_cuts_df.groupby(by='scaffold')

    no_cuts = intraexonic_cuts_count(exon_scaff_grouped, fp_cuts_grouped)
    print(
        f'After cut exon-breaking FPR is {100 * no_cuts / all_cuts_count:.2f}% as there are {all_cuts_count} total cuts'
    )

    print(f'--------------------------------Recall and exon-breaking precision ---------------------------------------')
    true_cuts = joined.query('cut == 1 and label == 1').shape[0]
    true_all = joined.query('label == 1')
    detectable_all = true_all.query('end - start <= 100 and end - start >= 40')

    print(f'Correctly cut {true_cuts}/{len(true_all)} ({strand} strand) introns.\n'
          f'Interfered with {no_cuts}/{exon_pos_df.shape[0]} ({strand} strand) exons.\n'
          f'Ratio {true_cuts / no_cuts:.2f}')

    recall = true_cuts / len(true_all) if true_cuts > 0 else 0
    exon_breaking_fpr = no_cuts / exon_pos_df.shape[0]

    print(f'Recall on detectable introns (40-100nt): {true_cuts}/{len(detectable_all)}')
    print(f'Total recall (all introns): {recall}')
    print(f'Exon breaking FPR: {exon_breaking_fpr}')

    return true_cuts, no_cuts, all_cuts_count, len(true_all), len(detectable_all), exon_pos_df.shape[0]


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
