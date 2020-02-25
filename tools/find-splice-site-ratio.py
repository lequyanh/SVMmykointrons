import sys

import pandas as pd
from Bio import SeqIO

from extract_tools import true_donorac_positions

"""
The purpose of this file is to determine how many GT dimers are true donors (given the fungi species)
Example output:

Number of positive strand introns: 11068
False splice site count: 2000451 
True splice site count: 10869.
Ratio: 184.05106265525808

"""


def main():
    # Introns file :/home/anhvu/Desktop/mykointrons-data/new-sequences/Blade1/Blade1-introns.fasta
    true_intron_file = sys.argv[1]
    # All donor candidates: e.g. Tripe1_results/splice-site-donor-dataset.csv
    intron_suspects_file = sys.argv[2]
    site = sys.argv[3]

    intron_suspects_positions = pd.read_csv(intron_suspects_file, delimiter=';')

    find_true_introns_proportion(true_intron_file, intron_suspects_positions, site)


# TODO rewrite to use panda merge
def find_true_introns_proportion(true_intron_file, intron_suspects_positions, site):
    with open(true_intron_file, 'r') as f:
        pos_strand_introns = list(filter(
            lambda seq_rec: seq_rec.description.split(' ')[1] == '+',
            SeqIO.parse(f, 'fasta')
        ))
        print(f'Number of positive strand introns: {len(pos_strand_introns)}')

        donor_positions, acceptor_positions = true_donorac_positions(true_intron_file)
        print('Number of true donors and acceptor with canonic GT/AG')
        print([sum(list(map(len, donor_positions.values()))), sum(list(map(len, acceptor_positions.values())))])

        false_splice_site_count = 0
        true_splice_site_count = 0

        exit(0)
        splicesite_positions_dict = donor_positions if site == 'donor' else acceptor_positions

        for _, row in intron_suspects_positions.iterrows():
            scaffold = row['scaffold']
            position = row['position']

            if scaffold not in splicesite_positions_dict.keys():
                false_splice_site_count += 1
                continue

            if position in splicesite_positions_dict[scaffold]:
                true_splice_site_count += 1
            else:
                false_splice_site_count += 1

        print(f'False splice site count: {false_splice_site_count} \nTrue splice site count: {true_splice_site_count}.'
              f'\nRatio: {false_splice_site_count / true_splice_site_count:.2f}')


if __name__ == "__main__":
    main()
