import csv
import logging
import sys

from fastalib import read_fasta

logging.basicConfig(
    level=logging.INFO,
    filename='purge-introns.log',
    filemode='w'
)

fasta_to_purge = sys.argv[1]  # '/home/anhvu/Desktop/Desktop/S21_contigs.fasta'
intron_locs = 'intron-result.csv'


def load_as_dicts(_fasta_to_purge: str, _intron_locs: str):
    # Load FASTA with DNA to cleanse introns from.
    # Load as dictionary where keys are scaffold names
    with open(_fasta_to_purge, 'r') as f:
        _scaffolds_input = {desc: seq for desc, seq in read_fasta(f)}

    # Load intron location data
    _intron_coords = dict()
    with open(_intron_locs, 'r') as f:
        # Assemble a scaffold:List[intron locations] dictionary
        for i, row in enumerate(csv.reader(f, delimiter=';')):
            if i == 0:
                continue  # Skip header of the csv

            scaff = row[0]
            positions = _intron_coords.get(scaff, [])
            positions.append((int(row[1]), int(row[2])))

            _intron_coords[scaff] = positions

    return _scaffolds_input, _intron_coords


def find_overlaps(_intron_coords: dict):
    # Process the introns (find overlap positions)
    _overlaps_dict = dict()
    for scaffold, positions in _intron_coords.items():
        last_end, last_start = 0, 0

        positions_non_overlap = []
        positions_overlap = []
        correction = 0

        for i, (start, end) in enumerate(positions):

            if start < last_end:
                positions_overlap += [(last_start, last_end, start, end)]

                if (last_start, last_end) in positions_non_overlap:
                    positions_non_overlap.remove((last_start, last_end))
                else:
                    # correction for multi-overlap
                    correction += 1

                overlap_ratio = (last_end - start) / (last_end - last_start)
                logging.info(f'{scaffold}-{last_start}-{last_end}---{start}-{end}-{overlap_ratio}')
            else:
                positions_non_overlap += [(start, end)]

            last_end = end
            last_start = start

        assert len(positions) == len(positions_non_overlap) + 2 * len(positions_overlap) - correction
        _overlaps_dict[scaffold] = (positions_non_overlap, positions_overlap)

    return _overlaps_dict


scaffold_sequences, intron_coords = load_as_dicts(fasta_to_purge, intron_locs)
overlaps_dict = find_overlaps(intron_coords)

for scaffold, (non_overlap_introns, overlap_introns) in overlaps_dict.items():
    scaffold_dna = scaffold_sequences[scaffold]
    print(scaffold)
