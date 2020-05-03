import logging
import os
import sys

import pandas as pd
from Bio import SeqIO, SeqRecord

logging.basicConfig(
    level=logging.INFO,
    filename='label-introns.log',
    filemode='w'
)


def main():
    # CSV file with intron candidates to label
    intron_candidates_csv = sys.argv[1]
    # FASTA file with true introns
    true_introns_fasta = sys.argv[2]

    len_range_start = int(sys.argv[3])
    len_range_end = int(sys.argv[4])

    if not os.path.isfile(true_introns_fasta):
        logging.warning(f'Introns file {true_introns_fasta} cannot be found. '
                        f'Candidates in {intron_candidates_csv} will be deleted')

        os.remove(intron_candidates_csv)
        exit(0)

    introns = get_introns_from_strand(true_introns_fasta, strand='+')

    right_length_introns = list(filter(lambda seq: len_range_start < len(seq) < len_range_end, introns))
    logging.info(f'Number of introns in the range {len_range_start} - {len_range_end}: {len(right_length_introns)}')

    intron_candidates = pd.read_csv(intron_candidates_csv, sep=';')

    label = [candidate['sequence'] in introns for i, candidate in intron_candidates.iterrows()]
    no_positive = sum(label)

    label = map(lambda flag: 1 if flag else -1, label)
    labeled_candidates = intron_candidates.assign(label=list(label))

    labeled_candidates.to_csv(intron_candidates_csv, sep=';', index=False)

    logging.info(f'{no_positive}/{len(intron_candidates)} candidates are positive introns')


def get_introns_from_strand(true_introns_fasta: str, strand: str):
    def is_strand(seq_rec: SeqRecord, sign: str):
        return seq_rec.description.split(' ')[1] == sign

    with open(true_introns_fasta, 'r') as intron_fasta_f:
        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))
        logging.info(f'Parsing introns from {true_introns_fasta} file.')
        logging.info(f'Total number of introns: {len(introns_seqrecords)}')

        introns = filter(lambda i_rec: is_strand(i_rec, strand), introns_seqrecords)
        introns = set([str(i.seq) for i in introns])

        logging.info(f'Number of true introns on the given (+) strand: {len(introns)}')

        return introns


if __name__ == '__main__':
    main()
