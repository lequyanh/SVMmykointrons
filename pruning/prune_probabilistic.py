import logging
import sys
from typing import List
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prune_tools import find_overlaps
from prune_tools import load_as_dicts
from prune_tools import load_length_counts_as_pdf
from prune_tools import prune_non_overlap_introns


def convert_coords(begin, end, start):
    intron = scaffold_sequences[scaffold][begin - 1:end]

    begin_converted = scaffold_prepruned.find(intron, start)
    end_converted = begin_converted + len(intron)

    return begin_converted, end_converted


def pdf_length_compare(begin1, end1, begin2, end2):
    length_1 = end1 - begin1
    length_2 = end2 - begin2
    p_len1, p_len2 = LENGTH_PDF.evaluate([length_1, length_2])

    return (begin1, end1) if p_len1 >= p_len2 else (begin2, end2)


def prune(
        scaffold_prepruned,
        _overlap_introns: List[Tuple[int, int, int, int]],
) -> List[str]:
    prev_fragment_origin, fragment_origin = 0, 0
    length_check = 0
    purged_fragments = []

    cut_coords, cut_length, pre_exon_length = None, None, None
    multioverlap_flag = False
    for begin1, end1, begin2, end2 in _overlap_introns:

        next_cut_coords = pdf_length_compare(begin1, end1, begin2, end2)

        if not cut_coords:
            begin_converted, old_begin = 0, 0
            cut_coords = next_cut_coords
            continue

        if multioverlap_flag:
            if cut_coords[1] > next_cut_coords[0]:
                print("Quadra overlap")
                fragment_origin = backup_fragment_origin
                length_check -= (cut_length + pre_exon_length)
            else:
                cut_coords = next_cut_coords  # Only for non-quadra overlap will be skipped
                multioverlap_flag = False
                continue

        multioverlap_flag = cut_coords[1] >= next_cut_coords[0]
        if multioverlap_flag:
            # re-evaluate the best cut
            cut_coords = pdf_length_compare(cut_coords[0], cut_coords[1], next_cut_coords[0], next_cut_coords[1])

        # Shift the area of search for intron as find() function finds only first occurrence.
        # In case of multi-overlap, when 3 introns can start at the same position, keep the old start
        find_start = begin_converted + 1 if old_begin < cut_coords[0] else begin_converted
        begin_converted, end_converted = convert_coords(cut_coords[0], cut_coords[1], find_start)

        if begin_converted == -1:
            logging.warning(f'Could not locate intron of coords {begin1, end1} or {begin2, end2} in '
                            f'pre-pruned scaffold {scaffold}')
            continue

        pre_exon_frags = scaffold_prepruned[fragment_origin:begin_converted]
        purged_fragments += [pre_exon_frags]

        cut_length = end_converted - begin_converted
        pre_exon_length = len(pre_exon_frags)
        length_check += (pre_exon_length + cut_length)

        backup_fragment_origin = fragment_origin
        fragment_origin = end_converted

        assert scaffold_prepruned[begin_converted:begin_converted + 2] == 'GT'
        assert scaffold_prepruned[end_converted - 2:end_converted] == 'AG'
        assert length_check == fragment_origin

        old_begin = cut_coords[0]
        cut_coords = next_cut_coords

    if _overlap_introns and not multioverlap_flag:
        find_start = begin_converted + 1 if old_begin < cut_coords[0] else begin_converted
        begin_converted, end_converted = convert_coords(cut_coords[0], cut_coords[1], find_start)

        pre_exon_frags = scaffold_prepruned[fragment_origin:begin_converted]
        purged_fragments += [pre_exon_frags]

        cut_length = end_converted - begin_converted
        pre_exon_length = len(pre_exon_frags)
        length_check += (pre_exon_length + cut_length)

        fragment_origin = end_converted

        assert length_check == fragment_origin

    if _overlap_introns:
        rest_of_exon = scaffold_prepruned[fragment_origin:]
        purged_fragments.append(rest_of_exon)

        length_check += len(rest_of_exon)
        assert length_check == len(scaffold_prepruned)

        return purged_fragments
    else:
        return [scaffold_prepruned]


if __name__ == "__main__":
    fasta_to_purge = sys.argv[1]  # /home/anhvu/Desktop/Desktop/S21_contigs.fasta
    intron_locs = sys.argv[2]  # CSV with intron locations
    counts_file = sys.argv[3]  # file with typical intron counts. Used to create prob. dens. function

    # fasta_to_purge = 'test'  # '/home/anhvu/Desktop/Desktop/S21_contigs.fasta'
    # intron_locs = 'positions.csv'
    # counts_file = 'sample.txt'

    scaffold_sequences, intron_coords = load_as_dicts(fasta_to_purge, intron_locs)
    LENGTH_PDF = load_length_counts_as_pdf(counts_file)
    overlaps_dict = find_overlaps(intron_coords)

    pruned_scaffolds = dict()
    for scaffold, (non_overlap_introns, overlap_introns) in overlaps_dict.items():
        print(f'Processing scaffold {scaffold}')
        scaffold_dna = scaffold_sequences[scaffold]

        scaffold_prepruned, exon_coord_mapping = prune_non_overlap_introns(scaffold_dna, non_overlap_introns)
        exon_fragments = prune(scaffold_prepruned, overlap_introns)

        # Each scaffold will be fragmented to exon pieces and overlap alternative windows
        pruned_scaffolds[scaffold] = exon_fragments

    with open(f'pruned-{fasta_to_purge}', 'w') as f:
        to_save = [SeqRecord(id=scaffold, seq=Seq(''.join(sequences))) for scaffold, sequences in
                   pruned_scaffolds.items()]
        SeqIO.write(to_save, f, 'fasta')
