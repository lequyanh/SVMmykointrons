import logging

from Bio import SeqIO

import fastalib as fl
import sys

# character that separates resulting columns
SEPARATOR = ';'

# assembly fasta file
assembly = sys.argv[1]
# donor dimer
donor = sys.argv[2]
# acceptor dimer
acceptor = sys.argv[3]
# window size for donors without the dimer: {lwindow}--donor--{rwindow}
donor_lwindow, donor_rwindow = int(sys.argv[4]), int(sys.argv[5])
# window size for acceptors without the dimer: {lwindow}--acceptor--{rwindow}
acceptor_lwindow, acceptor_rwindow = int(sys.argv[6]), int(sys.argv[7])
# strand
strand = sys.argv[8]

assert len(donor) == 2, 'donor must be a dimer'
assert len(acceptor) == 2, 'acceptor must be a dimer'
assert strand == 'plus' or strand == 'minus'

with open(assembly, 'r') as assembly_f:
    # loop through the whole assembly
    for scaffold_sr in SeqIO.parse(assembly_f, 'fasta'):
        scaffold = scaffold_sr.id
        sequence = str(scaffold_sr.seq)

        if strand == 'minus':
            sequence = fl.complementary(sequence)[::-1]  # Reverse complement

        for position, dimer in fl.dimers(sequence):
            if dimer == donor or dimer == acceptor:
                lwindow = donor_lwindow if dimer == donor else acceptor_lwindow
                rwindow = donor_rwindow if dimer == donor else acceptor_rwindow

                # compute start and end for the window
                win_start = position - lwindow - 1
                win_end = position + rwindow + 1

                if win_start < 0:
                    # Insufficient bases for window. Skipping
                    continue

                # extract the window
                window = sequence[win_start:win_end]

                # print scaffold, position and the window separated with the given separator
                print(scaffold, position, window, sep=SEPARATOR)
