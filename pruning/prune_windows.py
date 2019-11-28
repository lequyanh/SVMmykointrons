import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prune_tools import *

logging.basicConfig(
    level=logging.INFO,
    filename='prune-introns.log',
    filemode='w'
)


def convert_coords(begin, end, start):
    intron = scaffold_sequences[scaffold][begin - 1:end]

    begin_converted = scaffold_prepruned.find(intron, start)
    end_converted = begin_converted + len(intron)

    return begin_converted, end_converted


def prune(
        scaffold_prepruned,
        _overlap_introns: List[Tuple[int, int, int, int]],
) -> List[str]:

    prev_fragment_origin, fragment_origin = 0, 0
    length_check = 0
    purged_fragments = []
    begin1_converted, next_begin1_converted = None, None

    for begin1, end1, begin2, end2 in _overlap_introns:

        # Initialization for the first iteration
        if not begin1_converted:
            begin1_converted, end1_converted = convert_coords(begin1, end1, 0)
            begin2_converted, end2_converted = convert_coords(begin2, end2, 0)
            old_begin1 = begin1
            continue

        # Shift the area of search for intron as find() function finds only first occurrence.
        # In case of multi-overlap, when 3 introns can start at the same position, keep the old start
        find_start = begin1_converted + 1 if old_begin1 < begin1 else begin1_converted
        next_begin1_converted, next_end1_converted = convert_coords(begin1, end1, find_start)
        next_begin2_converted, next_end2_converted = convert_coords(begin2, end2, find_start)

        if next_begin1_converted == -1 or next_begin2_converted == -1:
            logging.warning(f'Could not locate intron of coords {begin1, end1} or {begin2, end2} in '
                            f'pre-pruned scaffold {scaffold}')
            continue

        max_end = max(end1_converted, end2_converted)

        multioverlap_flag = max_end > next_begin1_converted

        versions_start = max(fragment_origin, begin1_converted - WINDOW)
        if multioverlap_flag:
            # In case of multi-overlap we can't guarantee the right window won't affect the next overlap
            # We therefore set the conservative version end
            versions_end = max_end
        else:
            versions_end = min(max_end + WINDOW, next_begin1_converted)

        if versions_end - versions_start < 2 * WINDOW:
            versions_start = fragment_origin
            unaffected = ''
        else:
            # @unaffected is the part left from the current overlap (RNA from the last overlap to the current one)
            unaffected = scaffold_prepruned[fragment_origin:versions_start]

        pre_exon1 = scaffold_prepruned[versions_start:begin1_converted]
        pre_exon2 = scaffold_prepruned[versions_start:begin2_converted]

        post_exon1 = scaffold_prepruned[end1_converted:versions_end]
        post_exon2 = scaffold_prepruned[end2_converted:versions_end]

        if len(unaffected) < 2 * WINDOW:
            pre_exon1 = unaffected + pre_exon1
            pre_exon2 = unaffected + pre_exon2
            unaffected = ''
            purged_fragments += [pre_exon1 + post_exon1, pre_exon2 + post_exon2]
        else:
            purged_fragments += [unaffected, pre_exon1 + post_exon1, pre_exon2 + post_exon2]

        # Check for valid cut
        assert scaffold_prepruned[begin1_converted:begin1_converted + 2] == 'GT'
        assert scaffold_prepruned[end1_converted - 2:end1_converted] == 'AG'
        s1 = fragment_origin + len(unaffected) + len(pre_exon1) + (end1_converted - begin1_converted) + len(
            post_exon1)
        s2 = fragment_origin + len(unaffected) + len(pre_exon2) + (end2_converted - begin2_converted) + len(
            post_exon2)
        assert s1 == s2

        # check for multi-overlap
        # In case of 2-overlap the next @begin_converted should be way further, than the beginning of fragment.
        # If not, then multi-overlap occurred and we need to use the previous fragment origin
        if multioverlap_flag:
            assert s1 == max_end
            fragment_origin = prev_fragment_origin
        else:
            assert s1 == min(end2_converted + max(0, end1_converted - end2_converted) + WINDOW, next_begin1_converted)
            # Shift the @fragment_origin to the further of the two overlapped introns
            fragment_origin = max_end

        prev_fragment_origin = fragment_origin

        begin1_converted, end1_converted = next_begin1_converted, next_end1_converted
        begin2_converted, end2_converted = next_begin2_converted, next_end2_converted
        old_begin1 = begin1

    # Process the last overlap (only if there were any overlaps at all)
    if _overlap_introns:
        window_start = max(fragment_origin, begin1_converted - WINDOW)
        pre_exon1 = scaffold_prepruned[window_start:begin1_converted]
        pre_exon2 = scaffold_prepruned[window_start:begin2_converted]

        post_exon1 = scaffold_prepruned[end1_converted:end1_converted + WINDOW]
        post_exon2 = scaffold_prepruned[end2_converted:end2_converted + WINDOW]

        unaffected = scaffold_prepruned[fragment_origin:max(0, begin1_converted - WINDOW)]
        purged_fragments += [unaffected, pre_exon1 + post_exon1, pre_exon2 + post_exon2]

        # Now add the remainder of the sequence
        fragment_origin = max(end1_converted, end2_converted)

        rest_of_exon = scaffold_prepruned[fragment_origin:]
        purged_fragments.append(rest_of_exon)

        length_check += len(rest_of_exon)

        return purged_fragments
    else:
        return [scaffold_prepruned]


if __name__ == "__main__":
    fasta_to_purge = sys.argv[1]  # /home/anhvu/Desktop/Desktop/S21_contigs.fasta
    intron_locs = sys.argv[2]  # CSV with intron locations
    WINDOW = int(sys.argv[3])  # window around an overlap intron

    # fasta_to_purge = 'test'  # '/home/anhvu/Desktop/Desktop/S21_contigs.fasta'
    # intron_locs = 'positions.csv'
    # WINDOW = 10

    # Load FASTAS and intron positions as dictionaries, where the key is always the scaffold ID
    scaffold_sequences, intron_coords = load_as_dicts(fasta_to_purge, intron_locs)
    overlaps_dict = find_overlaps(intron_coords)

    pruned_scaffolds = dict()
    for scaffold, (non_overlap_introns, overlap_introns) in overlaps_dict.items():
        print(f'Processing scaffold {scaffold}')
        scaffold_dna = scaffold_sequences[scaffold]

        scaffold_prepruned, _ = prune_non_overlap_introns(scaffold_dna, non_overlap_introns)
        exon_fragments = prune(scaffold_prepruned, overlap_introns)

        # Each scaffold will be fragmented to exon pieces and overlap alternative windows
        pruned_scaffolds[scaffold] = exon_fragments

    with open(f'pruned-windows-{fasta_to_purge}', 'w') as f:
        to_save = []
        for scaffold, sequences in pruned_scaffolds.items():
            to_save += [SeqRecord(id=f'{scaffold}-frag{i}', seq=Seq(sequence)) for i, sequence in enumerate(sequences)]
        SeqIO.write(to_save, f, 'fasta')
