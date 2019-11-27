import csv
import logging
from typing import List
from typing import Tuple

from scipy import stats

from fastalib import read_fasta

logging.basicConfig(
    level=logging.INFO,
    filename='prune-introns.log',
    filemode='w'
)


def load_length_counts_as_pdf(counts_file: str):
    with open(counts_file, 'r') as f:
        counts = list(map(int, f.readlines()))
        kde = stats.gaussian_kde(counts)
    return kde


def load_as_dicts(fasta_to_purge: str, intron_locs: str):
    # Load FASTA with DNA to cleanse introns from.
    # Load as dictionary where keys are scaffold names
    with open(fasta_to_purge, 'r') as f:
        _scaffolds_input = {desc: seq for desc, seq in read_fasta(f)}

    # Load intron location data
    intron_coords = dict()
    with open(intron_locs, 'r') as f:
        # Assemble a scaffold:List[intron locations] dictionary
        for i, row in enumerate(csv.reader(f, delimiter=';')):
            if i == 0:
                continue  # Skip header of the csv

            scaff = row[0]
            positions = intron_coords.get(scaff, [])
            positions.append((int(row[1]), int(row[2])))

            intron_coords[scaff] = positions

    return _scaffolds_input, intron_coords


def find_overlaps(_intron_coords: dict):
    """
    Process the introns (find overlap positions)
    :return: Dictionary, where for each scaffold there are two lists - list of @non-overlap intron positions and
    list of @overlap_intron positions.
    """
    overlaps_dict = dict()
    for scaffold, positions in _intron_coords.items():
        last_end, last_start = 0, 0
        positions_non_overlap, positions_overlap = [], []
        correction = 0

        for i, (start, end) in enumerate(positions):

            if start < last_end:
                positions_overlap.append((last_start, last_end, start, end))

                if (last_start, last_end) in positions_non_overlap:
                    positions_non_overlap.remove((last_start, last_end))
                else:
                    correction += 1  # correction for multi-overlap

                overlap_ratio = (last_end - start) / (last_end - last_start)
                logging.info(f'{scaffold}-{last_start}-{last_end}---{start}-{end}-{overlap_ratio}')
            else:
                positions_non_overlap += [(start, end)]

            last_end = end
            last_start = start

        assert len(positions) == len(positions_non_overlap) + 2 * len(positions_overlap) - correction
        overlaps_dict[scaffold] = (positions_non_overlap, positions_overlap)

    return overlaps_dict


def prune_non_overlap_introns(
        scaffold_dna: str,
        non_overlap_introns: List[Tuple[int, int]]
) -> (str, List[Tuple[int, int]]):
    """
    Prunes a given DNA from introns, whose positions (start, end) are given in the @_non_overlap_introns list
    :return: Pruned scaffold and mapping between exon coordinates from the unpruned DNA to the pruned one
    """
    exon_begin = 0
    exon_coord_mapping = []
    purged_scaffold = ''
    for intron_begin, intron_end in non_overlap_introns:

        if intron_end - intron_begin > 80:
            continue  # Skipping long introns

        exon_end = intron_begin - 1  # Exon ends where the current intron begins. Starts where the previous intron ended
        exon_seq = scaffold_dna[exon_begin:exon_end]

        # Save coordinates. @exon_begin are the original ones. Coords in pruned DNA are just the current length
        # since the exons there are just glued together one by one.
        exon_coord_mapping.append((exon_begin, len(purged_scaffold)))
        purged_scaffold += exon_seq

        exon_begin = intron_end  # Shift the exon beginning to the intron end

    # Add the remainder of the sequence as exon (only if there has been any non-overlap introns)
    exon_coord_mapping.append((exon_begin, len(purged_scaffold)))
    exon_seq = scaffold_dna[exon_begin:]
    purged_scaffold += exon_seq

    return purged_scaffold, exon_coord_mapping


def convert_coords(
        exon_coord_mapping: List[Tuple[int, int]],
        overlap_intron_start: int,
        overlap_intron_end: int
):
    """"
    Converts the overlapped intron coordinates in the original FASTA to coordinates in pruned scaffolds.
    The overlapped introns lie in some of the exons.
    @:param _exon_coord_mapping: Map, with (exon_coordinates_unprunned : exon_coordinates_pruned)
    """

    def closest_smaller(cord_mapping):
        if overlap_intron_start > cord_mapping[0]:
            return overlap_intron_start - 1 - cord_mapping[0]
        return 9999999999999999999999

    # Find the coordinate mappings of the "bad" exon, that contain the overlapped intron
    bad_exon_start, bad_exon_converted_start = min(exon_coord_mapping, key=closest_smaller)

    begin_converted = bad_exon_converted_start + (overlap_intron_start - 1 - bad_exon_start)
    end_converted = bad_exon_converted_start + (overlap_intron_end - bad_exon_start)

    return begin_converted, end_converted
