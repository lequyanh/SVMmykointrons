import logging
import sys

import pandas as pd
from Bio import SeqIO

logging.getLogger().setLevel(logging.INFO)

if __name__ == '__main__':
    # CSV file with intron candidates to label
    intron_candidates_csv = sys.argv[1]
    # FASTA file with true introns
    true_introns_fasta = sys.argv[2]

    with open(true_introns_fasta, 'r') as intron_fasta_f:
        introns = SeqIO.parse(intron_fasta_f, 'fasta')
        introns = set([str(i.seq) for i in introns])

        i_candidates = pd.read_csv(intron_candidates_csv, sep=';')

        label = [row['sequence'] in introns for i, row in i_candidates.iterrows()]
        no_positive = sum(label)

        label = map(lambda flag: 1 if flag else -1, label)
        labeled_candidates = i_candidates.assign(label=list(label))

        labeled_candidates.to_csv(intron_candidates_csv, sep=';')

        logging.info(f'{no_positive}/{len(i_candidates)} candidates are positive introns')
