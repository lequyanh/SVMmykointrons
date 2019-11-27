import logging
import sys

import pandas as pd
from Bio import SeqIO, SeqRecord

logging.basicConfig(
    level=logging.INFO,
    filename='label-introns.log',
    filemode='w'
)

if __name__ == '__main__':
    # CSV file with intron candidates to label
    intron_candidates_csv = sys.argv[1]
    # FASTA file with true introns
    true_introns_fasta = sys.argv[2]

    with open(true_introns_fasta, 'r') as intron_fasta_f:
        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))
        logging.info(f'Total number of introns: {len(introns_seqrecords)}')

        # List which scaffolds the candidates come from
        intron_candidates = pd.read_csv(intron_candidates_csv, sep=';')
        candidate_scaffolds = intron_candidates.scaffold.unique()
        logging.info(f'Scaffolds, that the candidates were extracted from: {candidate_scaffolds}')

        # Keep only introns coming from scaffolds, that contain candidate introns
        def intron_scaffold(seq_rec: SeqRecord):
            return seq_rec.id.split(' ')[0]


        def is_strand(seq_rec: SeqRecord, sign: str):
            return seq_rec.description.split(' ')[1] == sign


        introns = filter(lambda i_rec: intron_scaffold(i_rec) in candidate_scaffolds and is_strand(i_rec, '+'),
                         introns_seqrecords)
        introns = set([str(i.seq) for i in introns])

        logging.info(f'Number of true introns on the given (+/-) strand: {len(introns)}')

        label = [candidate['sequence'] in introns for i, candidate in intron_candidates.iterrows()]
        no_positive = sum(label)

        label = map(lambda flag: 1 if flag else -1, label)
        labeled_candidates = intron_candidates.assign(label=list(label))

        labeled_candidates.to_csv(intron_candidates_csv, sep=';', index=False)

        logging.info(f'{no_positive}/{len(intron_candidates)} candidates are positive introns')
