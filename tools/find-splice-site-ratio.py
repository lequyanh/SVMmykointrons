import sys

import pandas as pd
from Bio import SeqIO

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
    # All donor candidates: e.g. Tripe1_results/splice-site-donor-results.csv
    intron_suspects_file = sys.argv[2]

    intron_suspects_positions = pd.read_csv(intron_suspects_file, delimiter=';')

    find_true_introns_proportion(true_intron_file, intron_suspects_positions)


def find_true_introns_proportion(true_intron_file, intron_suspects_positions):
    with open(true_intron_file, 'r') as f:
        pos_strand_introns = list(filter(lambda sr: sr.description.split(' ')[1] == '+', SeqIO.parse(f, 'fasta')))
        print(f'Number of positive strand introns: {len(pos_strand_introns)}')

        false_splice_site_count = 0
        true_splice_site_count = 0

        scaf_dict = dict()
        for sr in pos_strand_introns:
            descr_elements = sr.description.split(' ')
            scaff_idx = sr.id.split('_')[1]

            pos_arr = scaf_dict.get(scaff_idx, [])
            pos_arr.append(int(descr_elements[2]))
            scaf_dict[scaff_idx] = pos_arr

        for index, row in intron_suspects_positions.iterrows():
            scaffold_idx = row['scaffold'].split('_')[1]
            position = row['position']

            if scaffold_idx not in scaf_dict.keys():
                false_splice_site_count += 1
                continue

            if position in scaf_dict[scaffold_idx]:
                true_splice_site_count += 1
            else:
                false_splice_site_count += 1

        print(f'False splice site count: {false_splice_site_count} \nTrue splice site count: {true_splice_site_count}.'
              f'\nRatio: {false_splice_site_count / true_splice_site_count}')


if __name__ == "__main__":
    main()
