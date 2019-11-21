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
# intron_locs = 'intron-result.csv'
intron_locs = 'positions.csv'


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
    """
    Process the introns (find overlap positions)
    :param _intron_coords:
    :return: Dictionary, where for each scaffold there are two lists - list of non-overlap intron positions and
    list of overlap intron positions.
    """
    _overlaps_dict = dict()
    for scaffold, positions in _intron_coords.items():
        last_end, last_start = 0, 0
        positions_non_overlap, positions_overlap = [], []
        correction = 0

        for i, (start, end) in enumerate(positions):

            if start < last_end:
                positions_overlap.append([(last_start, last_end, start, end)])

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
        _overlaps_dict[scaffold] = (positions_non_overlap, positions_overlap)

    return _overlaps_dict


scaffold_sequences, intron_coords = load_as_dicts(fasta_to_purge, intron_locs)
overlaps_dict = find_overlaps(intron_coords)

purged_scaffolds = dict()
for scaffold, (non_overlap_introns, overlap_introns) in overlaps_dict.items():
    scaffold_dna = scaffold_sequences[scaffold]

    exon_begin = 0
    position_converter = []
    purged_scaffold = ''
    for intron_begin, intron_end in non_overlap_introns:
        exon_end = intron_begin - 1  # Exon ends where intron begins
        exon_seq = scaffold_dna[exon_begin:exon_end]

        position_converter.append([(exon_begin, len(purged_scaffold))])
        purged_scaffold += exon_seq

        exon_begin = intron_end  # Shift the exon beginning to the intron end

    # Add the remainder of the sequence as exon
    position_converter.append([(exon_begin, len(purged_scaffold))])
    exon_seq = scaffold_dna[exon_begin:]
    purged_scaffold += exon_seq

    purged_scaffolds[scaffold] = (purged_scaffold, position_converter)

    WINDOW = 100
    fragment_origin = 0
    for begin1, end1, begin2, end2 in overlap_introns:

        # version1
        bad_exon_start, bad_exon_converted_start = \
            min(position_converter, key=lambda original, purged: begin1 - original)

        begin1_converted = bad_exon_converted_start + (begin1 - bad_exon_start)  # start position of the
        pre_exon = purged_scaffold[begin1_converted - WINDOW:begin1_converted]

        end1_converted = bad_exon_converted_start + (end1 - bad_exon_start)
        post_exon = purged_scaffold[end1_converted:end1_converted + WINDOW]


        unaffected = purged_scaffold[fragment_origin:begin1_converted]
        # version2
        fragment_origin = end2_converted

