import csv
import itertools
from bisect import bisect

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.core.groupby import GroupBy

import fastalib as fl
import gfflib
from tools import performance_metrics


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
    shroom_name = 'Kocim1'
    gff_file = f'/home/anhvu/Desktop/mykointrons-data/new-sequences/{shroom_name}/Kocim1_GeneCatalog_genes_20160906.gff'
    assembly_fasta = f'/home/anhvu/Desktop/mykointrons-data/data/Assembly/{shroom_name}_AssemblyScaffolds.fasta'

    intron_annotations_file = f'intron-result.csv'
    cut_coords_file = f'cut_coords.csv'

    # Get the accuracy and recall of intron detection after the pruning step
    joined_df = join_on_position(cut_coords_file, intron_annotations_file)
    get_post_cut_accuracy_metrics(joined_df)

    # Examine, where false positive cuts fall into (exon or intergenic?)
    OUT_EXON_POS = f'{shroom_name}_exon_positions.csv'
    # First find where the exons are (their coordinates)
    extract_exon_positions(gff_file, assembly_fasta, OUT_EXON_POS)

    exon_pos_df = pd.read_csv(OUT_EXON_POS, delimiter=';')
    exon_scaff_grouped = exon_pos_df.groupby(by='scaffold')

    # See, where potential mistakes of pruning can be (intron dataset before classification)
    negative_intron_candidates = joined_df[joined_df.label == -1]
    intron_candidates_grouped = negative_intron_candidates.groupby(by='scaffold')
    print(f'--------------------- All false potential intra-exon cuts ---------------------------')
    no_cuts = determine_intraexonic_cuts(exon_scaff_grouped, intron_candidates_grouped)
    print(
        f'>> Cuts: {no_cuts} out of {len(negative_intron_candidates)} ({no_cuts / len(negative_intron_candidates):.2f})\n')

    # See, where false intron cuts happened
    false_positive_cuts_df = joined_df.query('pred == 1 and label == -1 and cut == 1')
    fp_cuts_grouped = false_positive_cuts_df.groupby(by='scaffold')

    print(f'--------------------- False positive intra-exon cuts ---------------------------')
    no_cuts = determine_intraexonic_cuts(exon_scaff_grouped, fp_cuts_grouped)
    print(f'>> Cuts: {no_cuts} out of {len(false_positive_cuts_df)} ({no_cuts / len(false_positive_cuts_df):.2f})')


def determine_intraexonic_cuts(exon_grouped: GroupBy, cuts_grouped: GroupBy) -> int:
    """
    Determine, how many (potential) cuts happened inside exon sequences.
    :param exon_grouped: GroupBy exon position for each scaffold
    :param cuts_grouped: GroupBy cut positions for each scaffold
    :return: Number of (potential) cuts inside exons
    """
    total_intraexon_cuts = 0
    for scaffold, group in cuts_grouped:
        try:
            scaff_exon_positions = exon_grouped.get_group(scaffold)
        except KeyError:
            print(f'No exons in scaffold {scaffold}, or scaffold missing')
            continue

        start_end_tuplelist = zip(scaff_exon_positions.start, scaff_exon_positions.end)
        start_end_tuplelist = sorted(start_end_tuplelist, key=lambda tup: tup[0])

        start_end_flattened = list(itertools.chain(*start_end_tuplelist))
        exon_length = sum([se[1] - se[0] for se in start_end_tuplelist])

        def inside_exon_cut(cut_start) -> bool:
            i = bisect(start_end_flattened, cut_start)
            return i % 2 == 1

        cut_starts = group['start']
        exon_cuts = [inside_exon_cut(cut_start) for cut_start in cut_starts]
        exon_cuts_no = sum(exon_cuts)
        print(f'{scaffold}: {exon_cuts_no} within exon cuts. '
              f'This makes {exon_cuts_no * 1000 / exon_length:.2f} cuts per 1k exon bp')

        total_intraexon_cuts += exon_cuts_no

    return total_intraexon_cuts


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


def extract_exon_positions(gff_file, assembly_fasta, output_csv: str) -> None:
    """
    Extracts exon coordinates into a csv with columns [scaffold, start, end]
    Only from positive strand
    :param gff_file: GFF file that contains fungi gene annotations (exons and their positions)
    :param assembly_fasta: FASTA file with whole fungi genome. Serves only as validation of extraction
    :param output_csv: CSV file with exon positions (scaffold and location)
    :return: None
    """
    db = gfflib.parse_gff(gff_file)

    with open(assembly_fasta, 'r') as f:
        scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    with open(output_csv, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=';')
        writer.writerow(['scaffold', 'start', 'end'])

        for gene in db.all_features(featuretype='gene'):
            for f in db.children(gene['name'][0], featuretype='exon'):
                # Extract the exon described in the GFF file from the given FASTA.
                exon_seq = f.sequence(assembly_fasta)

                scaffold = f.chrom
                strand = f.strand

                if strand == "+":
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert scaffold_seq[f.start - 1: f.end] == exon_seq

                    writer.writerow([scaffold, f.start - 1, f.end])


if __name__ == "__main__":
    main()
