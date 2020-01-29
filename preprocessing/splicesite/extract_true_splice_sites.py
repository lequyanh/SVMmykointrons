import logging
import sys

from Bio import SeqIO

from extract_tools import true_donorac_positions, extract_window

logging.basicConfig(
    level=logging.INFO,
    filename='extract-true-windows.log',
    filemode='w'
)


def main():
    shroom_name = sys.argv[1]
    assembly_folder = sys.argv[2]
    introns_folder = sys.argv[3]

    margin_size = int(sys.argv[4])
    # ---------------------------
    # shroom_name = 'Ramac1'
    # assembly_folder = "/home/anhvu/Desktop/mykointrons-data/data/Assembly"
    # introns_folder = "/home/anhvu/Desktop/mykointrons-data/new-sequences"
    #
    # margin_size = 200

    assembly_fasta = f'{assembly_folder}/{shroom_name}_AssemblyScaffolds.fasta'
    introns_fasta = f'{introns_folder}/{shroom_name}/{shroom_name}-introns.fasta'

    donor_positions, acceptor_positions = true_donorac_positions(introns_fasta)

    print(f'Writing true donor/acceptor windows (on positive strand) for {shroom_name}')

    true_donors, true_acceptors = \
        retrieve_true_splice_sites(assembly_fasta, donor_positions, acceptor_positions, margin_size)

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-donor-true.fasta', 'w') as f:
        SeqIO.write(true_donors, f, 'fasta')

    with open(f'{introns_folder}/{shroom_name}/{shroom_name}-acceptor-true.fasta', 'w') as f:
        SeqIO.write(true_acceptors, f, 'fasta')


def retrieve_true_splice_sites(
        assembly_fasta: str,
        donor_positions: dict,
        acceptor_positions: dict,
        margin_size: int
):
    true_donors = list()
    true_acceptors = list()

    with open(assembly_fasta, 'r') as assembly_f:
        # loop through the whole assembly
        for seq_rec in SeqIO.parse(assembly_f, 'fasta'):
            scaffold = seq_rec.description
            sequence = str(seq_rec.seq)

            for donor_pos in donor_positions.get(scaffold, []):
                window = extract_window(sequence, donor_pos, margin_size, scaffold)
                true_donors = true_donors + [window] if window else true_donors

            for acc_pos in acceptor_positions.get(scaffold, []):
                window = extract_window(sequence, acc_pos, margin_size, scaffold)
                true_acceptors = true_acceptors + [window] if window else true_acceptors

        print(f'\t True splice sites: {len(true_donors)} donor, {len(true_acceptors)} acceptor windows')

        return true_donors, true_acceptors


if __name__ == "__main__":
    main()
