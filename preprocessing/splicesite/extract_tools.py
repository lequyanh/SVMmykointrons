import logging
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DONOR = 'GT'
ACCEPTOR = 'AG'


def extract_window(
        sequence: str,
        position: int,
        margin_size: int,
        scaffold: str
) -> Optional[SeqRecord]:
    # compute start and end for the window
    win_start = position - margin_size - 1
    win_end = position + margin_size + 1

    if win_start < 0:
        return None

    # extract the window
    window = sequence[win_start:win_end]

    dimer = window[margin_size:margin_size + 2]
    if dimer not in {DONOR, ACCEPTOR}:
        return None

    if len(window) < 2 * margin_size + 2:
        logging.warning("Windows too short, skipping...")
        return None

    return SeqRecord(
        id=' '.join([scaffold, '+', str(win_start), str(win_end)]),
        seq=Seq(window)
    )


def true_donorac_positions(introns_fasta: str):
    with open(introns_fasta, 'r') as intron_fasta_f:
        donor_positions = dict()
        acceptor_positions = dict()

        non_canonic = 0
        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))
        for isr in introns_seqrecords:
            scaff_id, strand, start, end = isr.description.split(' ')

            if strand == '+':
                donor_start = int(start)
                acceptor_start = int(end) - 1

                don_pos = donor_positions.get(scaff_id, [])
                acc_pos = acceptor_positions.get(scaff_id, [])

                don_pos.append(donor_start)
                acc_pos.append(acceptor_start)

                donor_positions[scaff_id] = don_pos
                acceptor_positions[scaff_id] = acc_pos

                if str(isr.seq[:2]) not in {DONOR, ACCEPTOR}:
                    non_canonic += 1

        print(f'Parsed intron file. '
              f'Non-canonic dimers: {non_canonic} out of 2x{len(introns_seqrecords)} (num of introns)')

    return donor_positions, acceptor_positions

#
# def get_introns_position_info(introns_locs_fasta: str) -> Iterable[List]:
#     """
#     Parses intron file. Returns list, where each element is a intron "object"
#     This object is a list of scaffold, position and sequence
#     :param introns_locs_fasta: Location of introns file
#     """
#     # Fasta file with intron locations and sequences
#     with open(introns_locs_fasta, 'r') as introns_locs_f:
#         lines = introns_locs_f.readlines()
#
#         # From 1st line every second line (scaffold, strand, start, end)
#         intron_descriptions = lines[::2]
#         # From 2nd line every second line (intron sequence)
#         intron_seqs = lines[1::2]
#
#         def join_sequence_with_info(rec):
#             intron_info = rec[0].rstrip().split(' ')
#             intron_seq = rec[1].rstrip()
#
#             intron_info.append(intron_seq)
#             return intron_info
#
#         # Join parsed location info with the sequence itself
#         intron_descriptions = map(
#             join_sequence_with_info,
#             zip(intron_descriptions, intron_seqs)
#         )
#
#         return intron_descriptions
