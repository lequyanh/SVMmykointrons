import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prune_tools import *


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
    purged_fragments, introns = [], []

    def cut_region(fragment_start: int, cut_idxs: [int, int]):
        begin_converted, end_converted = convert_coords(cut_idxs[0], cut_idxs[1], fragment_start)

        if begin_converted == -1:
            begin_converted, end_converted = convert_coords(cut_idxs[0], cut_idxs[1], fragment_start - 1)

        pre_intron_region = scaffold_prepruned[fragment_start:begin_converted]
        intron = scaffold_prepruned[begin_converted:end_converted]

        purged_fragments.append(pre_intron_region)
        introns.append(intron)

        length_check = len(''.join(purged_fragments)) + len(''.join(introns))

        prev_fragment_origin = fragment_start
        fragment_start = end_converted

        assert intron[:2] == 'GT' and intron[- 2:] == 'AG'
        assert length_check > len(''.join(purged_fragments))
        try:
            assert length_check == end_converted
        except AssertionError:
            logging.warning(f'Fragment origin is {fragment_start} but length check is {length_check} '
                            f'Scaffold {scaffold}')

        return purged_fragments, introns, fragment_start, prev_fragment_origin

    cut_coords = None
    multioverlap_flag = False
    for begin1, end1, begin2, end2 in _overlap_introns:

        next_cut_coords = pdf_length_compare(begin1, end1, begin2, end2)

        if not cut_coords:
            begin_converted, old_begin = 0, 0
            cut_coords = next_cut_coords
            continue

        if multioverlap_flag:
            if cut_coords[1] >= next_cut_coords[0]:
                logging.info("Triple overlap")
                fragment_origin = prev_fragment_origin
                # noinspection PyUnboundLocalVariable
                purged_fragments.pop()
                # noinspection PyUnboundLocalVariable
                introns.pop()
            else:
                cut_coords = next_cut_coords  # Only for non-quadra overlap will be skipped
                multioverlap_flag = False
                continue

        multioverlap_flag = cut_coords[1] >= next_cut_coords[0]
        # In case of double overlap, re-evaluate the best cut
        if multioverlap_flag:
            cut_coords = pdf_length_compare(cut_coords[0], cut_coords[1], next_cut_coords[0], next_cut_coords[1])

        print(f'{scaffold};{cut_coords[0]};{cut_coords[1]}')
        purged_fragments, introns, fragment_origin, prev_fragment_origin = cut_region(fragment_origin, cut_coords)
        cut_coords = next_cut_coords

    if _overlap_introns and not multioverlap_flag:
        # Process the last overlap (if not multi-overlap)
        purged_fragments, introns, fragment_origin, _ = cut_region(fragment_origin, cut_coords)

    if _overlap_introns:
        # Processed all overlaps. Just add the rest of the sequence
        rest_of_exon = scaffold_prepruned[fragment_origin:]
        purged_fragments.append(rest_of_exon)

        try:
            assert len(''.join(purged_fragments)) + len(''.join(introns)) == len(scaffold_prepruned)
        except AssertionError:
            logging.warning(
                f'Fragment length is {len(scaffold_prepruned)} but length check is {len("".join(purged_fragments)) + len("".join(introns))} '
                f'Scaffold {scaffold}')

        return purged_fragments
    else:
        return [scaffold_prepruned]


if __name__ == "__main__":
    fasta_to_purge = sys.argv[1]  # /home/anhvu/Desktop/Desktop/S21_contigs.fasta
    intron_locs = sys.argv[2]  # CSV with intron locations
    counts_file = sys.argv[3]  # file with typical intron counts. Used to create prob. dens. function
    # /home/anhvu/PycharmProjects/mycointrons/statistics/intron-lenghts-statistics/basidiomycota-intron-lens.txt

    scaffold_sequences, intron_coords = load_as_dicts(fasta_to_purge, intron_locs)
    LENGTH_PDF = load_length_counts_as_pdf(counts_file)
    overlaps_dict = find_overlaps(intron_coords)

    no_intron_scaffolds = get_scaffold_without_introns(scaffold_sequences, overlaps_dict)  # type: dict

    pruned_scaffolds = dict()
    for scaffold, (non_overlap_introns, overlap_introns) in overlaps_dict.items():
        logging.info(f'Processing scaffold {scaffold}')
        scaffold_dna = scaffold_sequences[scaffold]

        scaffold_prepruned, _ = prune_non_overlap_introns(scaffold, scaffold_dna, non_overlap_introns)
        exon_fragments = prune(scaffold_prepruned, overlap_introns)

        pruned_scaffolds[scaffold] = ''.join(exon_fragments)
        assert len(scaffold_dna) >= len(pruned_scaffolds[scaffold])

    assert set(pruned_scaffolds.keys()).intersection(set(no_intron_scaffolds.keys())) == set()
    assert scaffold_sequences.keys() == {**pruned_scaffolds, **no_intron_scaffolds}.keys()
    assert len(''.join(pruned_scaffolds.values())) + len(''.join(no_intron_scaffolds.values())) < len(
        ''.join(scaffold_sequences.values()))

    filename = os.path.basename(fasta_to_purge)
    with open(f'pruned-{filename}', 'w') as f:
        # Create SeqRecords from joined dictionaries - purged scaffolds and scaffolds with no introns
        to_save = [
            SeqRecord(id=scaffold, seq=Seq(sequence))
            for scaffold, sequence in
            {**pruned_scaffolds, **no_intron_scaffolds}.items()
        ]

        SeqIO.write(to_save, f, 'fasta')
