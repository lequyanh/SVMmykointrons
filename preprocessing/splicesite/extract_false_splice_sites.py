import logging
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import fastalib as fl

logging.basicConfig(
    level=logging.INFO,
    filename='extract-false-windows.log',
    filemode='w'
)

MARGIN_SIZE = 150
DONOR = 'GT'
ACCEPTOR = 'AG'


def retrieve_true_intron_positions(introns_fasta: str):
    with open(introns_fasta, 'r') as intron_fasta_f:
        donor_positions = dict()
        acceptor_positions = dict()

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

    return donor_positions, acceptor_positions


def retrieve_false_splice_sites(donor_positions: dict, acceptor_positions: dict):
    false_donors = list()
    false_acceptors = list()
    k = 0

    with open(assembly_fasta, 'r') as assembly_f:
        # loop through the whole assembly
        for scaffold, sequence in fl.read_fasta(assembly_f):
            for position, dimer in fl.dimers(sequence):
                if dimer == DONOR or dimer == ACCEPTOR:

                    # compute start and end for the window
                    win_start = position - MARGIN_SIZE - 1
                    win_end = position + MARGIN_SIZE + 1

                    if win_start < 0:
                        logging.info(f'Window is too short')
                        continue

                    if dimer == DONOR:
                        if position in donor_positions.get(scaffold, []):
                            logging.info(f'Dimer is a splice site, skip it')
                            continue
                    else:
                        if position in acceptor_positions.get(scaffold, []):
                            logging.info(f'Dimer is a splice site, skip it')
                            continue

                    # extract the window
                    window = sequence[win_start:win_end]
                    assert window[MARGIN_SIZE:MARGIN_SIZE + 2] in ['GT', 'AG']

                    if len(window) < 2 * MARGIN_SIZE + 2:
                        logging.info(f'Window too short')
                        continue

                    false_splice_site = SeqRecord(
                        id=' '.join([scaffold, '+', str(win_start), str(win_end)]),
                        seq=Seq(window)
                    )

                    if dimer == DONOR:
                        false_donors.append(false_splice_site)
                    else:
                        false_acceptors.append(false_splice_site)

                    k += 1
                    if k >= false_examples_limit:
                        logging.info(f'Splice sites written: {k}')
                        return false_donors, false_acceptors

        logging.info(f'Splice sites written: {k}')
        return false_donors, false_acceptors


if __name__ == '__main__':
    shroom_name = sys.argv[1]
    assembly_folder = sys.argv[2]
    introns_folder = sys.argv[3]
    false_examples_limit = int(sys.argv[4])
    # ---------------------------
    # shroom_name = 'Aciaci1'
    # assembly_folder = "/home/anhvu/Desktop/mykointrons-data/data/Assembly"
    # introns_folder = "/home/anhvu/Desktop/mykointrons-data/new-sequences"
    assembly_fasta = f'{assembly_folder}/{shroom_name}_AssemblyScaffolds.fasta'
    _introns_fasta = f'{introns_folder}/{shroom_name}/{shroom_name}-introns.fasta'

    _donor_positions, _acceptor_positions = retrieve_true_intron_positions(_introns_fasta)
    _false_donors, _false_acceptors = retrieve_false_splice_sites(_donor_positions, _acceptor_positions)

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-donor-false2.fasta', 'w') as f:
        SeqIO.write(_false_donors, f, 'fasta')

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-acceptor-false2.fasta', 'w') as f:
        SeqIO.write(_false_acceptors, f, 'fasta')
