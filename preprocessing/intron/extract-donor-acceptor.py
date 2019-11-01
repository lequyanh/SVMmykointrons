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
# number of scaffolds to consider
no_scaffolds = int(sys.argv[8])

assert len(donor) == 2, 'donor must be a dimer'
assert len(acceptor) == 2, 'acceptor must be a dimer'

with open(assembly, 'r') as assembly_f:
    # Load scaffolds in assembly
    scaffolds = list(SeqIO.parse(assembly_f, 'fasta'))
    logging.info(f'Loaded {len(scaffolds)} scaffolds in {assembly} FASTA file.')

    # Sample a set of assemblies
    sel_scaff_idx = np.random.choice(len(scaffolds), no_scaffolds)
    selected_scaffs = np.array(scaffolds)[sel_scaff_idx]

    for scaffold in selected_scaffs:
        sequence = scaffold.seq
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
                scaffold_info = f'{shroom}:{scaffold.id}'
                print(scaffold_info, position, window, sep=SEPARATOR)
