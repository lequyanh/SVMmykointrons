import logging
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO

import fastalib as fl

logging.getLogger().setLevel(logging.INFO)

# character that separates resulting columns
SEPARATOR = ';'

# assembly fasta file
assembly = sys.argv[1]
# donor dimer
donor = sys.argv[2]
# acceptor dimer
acceptor = sys.argv[3]
# window size for donors without the dimer: {   lwindow}--donor--{rwindow}
donor_lwindow, donor_rwindow = int(sys.argv[4]), int(sys.argv[5])
# window size for acceptors without the dimer: {lwindow}--acceptor--{rwindow}
acceptor_lwindow, acceptor_rwindow = int(sys.argv[6]), int(sys.argv[7])
# max number of scaffolds splice site candidates to consider
splice_sites_limit = int(sys.argv[8])
# strand
strand = sys.argv[9]

assert len(donor) == 2, 'donor must be a dimer'
assert len(acceptor) == 2, 'acceptor must be a dimer'
assert strand == '+' or strand == '-'

with open(assembly, 'r') as assembly_f:
    # Load scaffolds in assembly
    scaffolds = list(SeqIO.parse(assembly_f, 'fasta'))
    no_scaffolds = len(scaffolds)
    logging.info(f'Loaded {no_scaffolds} scaffolds in {assembly} FASTA file.')

    # Sample a set of assemblies
    scaffolds = np.array(scaffolds)

    np.random.seed(42)
    np.random.shuffle(scaffolds)

    k = 0

    for i, scaffold in enumerate(scaffolds):
        sequence = scaffold.seq if strand == '+' else scaffold.seq.reverse_complement()

        for position, dimer in fl.dimers(sequence):
            if dimer == donor or dimer == acceptor:
                lwindow = donor_lwindow if dimer == donor else acceptor_lwindow
                rwindow = donor_rwindow if dimer == donor else acceptor_rwindow

                # compute start and end for the window
                win_start = position - lwindow - 1
                win_end = position + rwindow + 1

                # extract the window
                window = sequence[win_start:win_end]

                # print scaffold, position and the window separated with the given separator
                shroom = Path(assembly).parts[-1][:-6]
                print(scaffold.id, position, window, sep=SEPARATOR)

                k += 1
                # We don't want a donor to be without a candidate acceptor - must terminate on acceptor
                # Can be enhanced by keeping track of last position, so we terminate on an acceptor in reasonable dist.
                if k >= splice_sites_limit and dimer == acceptor:
                    logging.info(f'Splice sites written {k} from {i} scaffolds')
                    exit(0)

    logging.info(f'Splice sites written {k} from {i} scaffolds')
