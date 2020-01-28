import logging
import sys

from Bio import SeqIO

import fastalib as fl
from extract_tools import true_donorac_positions, DONOR, ACCEPTOR, extract_window

logging.basicConfig(
    level=logging.INFO,
    filename='extract-false-windows.log',
    filemode='w'
)


def main():
    # shroom_name = sys.argv[1]
    # assembly_folder = sys.argv[2]
    # introns_folder = sys.argv[3]
    #
    # examples_limit = int(sys.argv[4])
    # margin_size = int(sys.argv[5])
    # ---------------------------
    shroom_name = 'Lenvul1'
    assembly_folder = "/home/anhvu/Desktop/mykointrons-data/data/Assembly"
    introns_folder = "/home/anhvu/Desktop/mykointrons-data/new-sequences"

    examples_limit = 50000
    margin_size = 200

    assembly_fasta = f'{assembly_folder}/{shroom_name}_AssemblyScaffolds.fasta'
    introns_fasta = f'{introns_folder}/{shroom_name}/{shroom_name}-introns.fasta'

    donor_positions, acceptor_positions = true_donorac_positions(introns_fasta)
    false_donors, false_acceptors = \
        retrieve_false_splice_sites(assembly_fasta, donor_positions, acceptor_positions, margin_size, examples_limit)

    print(f'writing false windows for {shroom_name}')

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-donor-false.fasta', 'w') as f:
        SeqIO.write(false_donors, f, 'fasta')

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-acceptor-false.fasta', 'w') as f:
        SeqIO.write(false_acceptors, f, 'fasta')


def retrieve_false_splice_sites(
        assembly_fasta: str,
        donor_positions: dict,
        acceptor_positions: dict,
        margin_size: int,
        examples_limit: int
):
    false_donors = list()
    false_acceptors = list()
    k = 0

    def yield_false_windows():
        print(f'False splice sites written: {k}/{examples_limit}. '
              f'{len(false_donors)} donor, {len(false_acceptors)} acceptor windows')
        return false_donors, false_acceptors

    with open(assembly_fasta, 'r') as assembly_f:

        for scaffold, sequence in fl.read_fasta(assembly_f):
            for position, dimer in fl.dimers(sequence):
                if k >= examples_limit:
                    return yield_false_windows()

                if dimer == DONOR and position not in donor_positions.get(scaffold, []):
                    window = extract_window(sequence, position, margin_size, scaffold)
                    false_donors = false_donors + [window] if window else false_donors
                    k += 1

                elif dimer == ACCEPTOR and position not in acceptor_positions.get(scaffold, []):
                    window = extract_window(sequence, position, margin_size, scaffold)
                    false_acceptors = false_acceptors + [window] if window else false_acceptors
                    k += 1
                else:
                    logging.info(f'Dimer is a splice site, skip it')
                    continue

        return yield_false_windows()


if __name__ == '__main__':
    main()
