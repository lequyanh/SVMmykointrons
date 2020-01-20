import pandas as pd
from Bio import SeqIO

splice_sites_windows = pd.read_csv('Tripe1_results/splice-site-donor-dataset.csv', delimiter=';')

with open("/home/anhvu/Desktop/mykointrons-data/new-sequences/Blade1/Blade1-introns.fasta", 'r') as f:
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

    for index, row in splice_sites_windows.iterrows():
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
