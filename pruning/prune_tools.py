import logging
from typing import List
from typing import Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    filename='prune-tools.log',
    filemode='w'
)


def load_length_counts_as_pdf(counts_file: str):
    with open(counts_file, 'r') as f:
        counts = list(map(int, f.readlines()))
        kde = stats.gaussian_kde(counts)
    return kde


def load_as_dicts(fasta_to_purge: str, intron_locs: str, strand: str):
    # Load FASTA with DNA to cleanse introns from.
    with open(fasta_to_purge, 'r') as f:
        sequences = SeqIO.parse(f, 'fasta')

        if strand == 'minus':
            sequences = [SeqRecord(sr.seq.reverse_complement(), sr.id) for sr in sequences]

        # Load as dictionary where keys are scaffold names
        scaffolds_input = {seq_record.id: str(seq_record.seq) for seq_record in sequences}

    # Load intron location data
    intron_pos_df = pd.read_csv(intron_locs, delimiter=";")
    if 'prediction' in intron_pos_df.columns:
        intron_pos_df_true = intron_pos_df[intron_pos_df['prediction'] == 1]
    else:
        intron_pos_df_true = intron_pos_df

    grouped_intron_pos = intron_pos_df_true.groupby(by='scaffold')
    intron_coords = dict()
    for scaffold, positions in grouped_intron_pos:
        intron_coords[scaffold] = positions.apply(lambda x: (x.start, x.end), axis=1).values

    return scaffolds_input, intron_coords


def find_overlaps(_intron_coords: dict):
    """
    Process the introns (find overlap positions)
    :return: Dictionary, where for each scaffold there are two lists - list of @non-overlap intron positions and
    list of @overlap_intron positions.
    """
    overlaps_dict = dict()
    for scaffold, positions in _intron_coords.items():
        # The "last last" vars are needed in case the overlap is not between two immediate neighbours (it skips one)
        # e.g. for intron positions 10-15, 12-14, 14-16)
        last_end, last_start, max_end_start, max_end = 0, 0, 0, 0
        positions_non_overlap, positions_overlap = [], []
        correction = 0

        for i, (start, end) in enumerate(positions):

            if start <= last_end:
                positions_overlap.append((last_start, last_end, start, end))

                if (last_start, last_end) in positions_non_overlap:
                    positions_non_overlap.remove((last_start, last_end))
                else:
                    correction += 1  # correction for multi-overlap

                overlap_ratio = (last_end - start) / (last_end - last_start)
                logging.info(f'{scaffold}> simple overlap ratio {overlap_ratio}')
            elif last_end <= start <= max_end:
                # remove the last overlap. That intron is too short anyways, so ignore it
                positions_overlap = positions_overlap[:-1]
                logging.warning("Overlap type: skip one")

                positions_overlap.append((max_end_start, max_end, start, end))
                # even though we should increment for multi-overlap count, we removed one overlap pair
                # Therefore we removed 2 intron positions and added one overlap position, i.e. -2 + 1
                correction -= 1
            else:
                positions_non_overlap += [(start, end)]

            max_end = max(max_end, end)
            last_end = end
            max_end_start = start if max_end == end else max_end_start
            last_start = start

        assert len(positions) == len(positions_non_overlap) + 2 * len(positions_overlap) - correction
        overlaps_dict[scaffold] = (positions_non_overlap, positions_overlap)

    return overlaps_dict


def prune_non_overlap_introns(
        scaffold: str,
        scaffold_dna: str,
        non_overlap_introns: List[Tuple[int, int]]
) -> (str, List[Tuple[int, int]]):
    """
    Prunes a given DNA from introns, whose positions (start, end) are given in the @_non_overlap_introns list
    :return: Pruned scaffold and mapping between exon coordinates from the unpruned DNA to the pruned one
    """
    exon_begin = 0
    purged_scaffold = ''
    for intron_begin, intron_end in non_overlap_introns:

        if intron_end - intron_begin > 80:
            continue  # Skipping long introns

        print(f'{scaffold};{intron_begin};{intron_end}')

        exon_end = intron_begin - 1  # Exon ends where the current intron begins. Starts where the previous intron ended
        exon_seq = scaffold_dna[exon_begin:exon_end]

        # Save coordinates. @exon_begin are the original ones. Coords in pruned DNA are just the current length
        # since the exons there are just glued together one by one.
        purged_scaffold += exon_seq

        exon_begin = intron_end  # Shift the exon beginning to the intron end

    # Add the remainder of the sequence as exon (only if there has been any non-overlap introns)
    exon_seq = scaffold_dna[exon_begin:]
    purged_scaffold += exon_seq

    return purged_scaffold


def get_scaffold_without_introns(original_scaffolds: dict, to_prune_scaffolds: dict):
    scaffold_no_intron = set(original_scaffolds.keys()) - set(to_prune_scaffolds.keys())

    return {sf: original_scaffolds[sf] for sf in scaffold_no_intron}
