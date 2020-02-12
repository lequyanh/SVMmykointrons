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
    @cut_coords_file  Prune introns and record coords which were cut. Can be obtained by running
                            echo "start;end" > cut_coords.csv &
                            python ../pruning/prune_probabilistic.py
                                    ~/Desktop/mykointrons-data/data/Assembly/Kocim1_AssemblyScaffolds.fasta
                                    intron-result.csv
                                    ../statistics/intron-lenghts-statistics/basidiomycota-intron-lens.txt
                            >> cut_coords.csv
    @intron_annotation_file     Run pipeline with an exception - do not discard negatively labeled introns candidates
    """
    shroom_name = 'Tripe1'
    gff_file = f'/home/anhvu/Desktop/mykointrons-data/new-sequences/{shroom_name}/Tripe1_GeneCatalog_genes_20150117.gff'
    assembly_fasta = f'/home/anhvu/Desktop/mykointrons-data/data/Assembly/{shroom_name}_AssemblyScaffolds.fasta'
    intron_annotations_file = f'{shroom_name}_intron_result.csv'
    cut_coords_file = f'{shroom_name}_cut_coords.csv'

    # Get the accuracy and recall after the pruning step
    joined_df = join_on_position(cut_coords_file, intron_annotations_file)
    get_post_cut_accuracy_metrics(joined_df)

    # Examine, where false positive cuts fall (exon or intergenic?)
    OUT_EXON_POS = f'{shroom_name}_exon_positions.csv'
    extract_exon_positions(gff_file, assembly_fasta, OUT_EXON_POS)

    exon_pos_df = pd.read_csv(OUT_EXON_POS, delimiter=';')
    exon_scaff_grouped = exon_pos_df.groupby(by='scaffold')

    negative_intron_candidates = joined_df[joined_df.label == -1]
    intron_candidates_grouped = negative_intron_candidates.groupby(by='scaffold')
    print(f'--------------------- All false potential intra-exon cuts ---------------------------')
    no_cuts = determine_intraexonic_cuts(exon_scaff_grouped, intron_candidates_grouped)
    print(f'>> Cuts: {no_cuts} out of {len(negative_intron_candidates)} ({no_cuts / len(negative_intron_candidates):.2f})\n')

    false_positive_cuts_df = joined_df.query('pred == 1 and label == -1 and cut == 1')
    fp_cuts_grouped = false_positive_cuts_df.groupby(by='scaffold')

    print(f'--------------------- False positive intra-exon cuts ---------------------------')
    no_cuts = determine_intraexonic_cuts(exon_scaff_grouped, fp_cuts_grouped)
    print(f'>> Cuts: {no_cuts} out of {len(false_positive_cuts_df)} ({no_cuts / len(false_positive_cuts_df):.2f})')


def determine_intraexonic_cuts(exon_grouped: GroupBy, cuts_grouped: GroupBy):
    total_interexon_cuts = 0
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

        total_interexon_cuts += exon_cuts_no

    return total_interexon_cuts


def get_post_cut_accuracy_metrics(joined):
    metrics_before = performance_metrics(joined['label'], joined['pred'], None)
    metrics_before = '\n'.join(metrics_before)

    print('Before cut metrics')
    print(metrics_before)
    print('\n')

    metrics = performance_metrics(joined['label'], joined['cut'], None)
    metrics = '\n'.join(metrics)

    print('After cut metrics')
    print(metrics)


def join_on_position(cut_coords_file, intron_annotations_file) -> DataFrame:
    cut_coords_df = pd.read_csv(cut_coords_file, delimiter=';')
    intron_coords_df = pd.read_csv(intron_annotations_file, delimiter=';')

    cut_coords_df['cut'] = np.ones(len(cut_coords_df))

    joined = intron_coords_df.merge(cut_coords_df, how="left")
    joined[['cut']] = joined[['cut']].fillna(value=-1)

    return joined


def extract_exon_positions(gff_file, assembly_fasta, output_csv: str):
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
