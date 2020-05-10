import sys

import numpy as np
from Bio import SeqIO

ASSEMBLIES_LOC = '/home/anhvu/Desktop/mykointrons-data/data/Assembly'
NEWSEQUENCES_LOC = '/home/anhvu/Desktop/mykointrons-data/new-sequences'


def main():
    # assembly fasta file
    fungi_name = sys.argv[1]
    # portion of scaffolds to take
    out_reduced = sys.argv[2]  # /home/anhvu/Desktop/mykointrons-data/reduced/
    scaffolds_portion = float(sys.argv[3])  # 0.04

    extract_fungi_portion(fungi_name, scaffolds_portion, out_reduced)


def main2():
    fungi_names_file = sys.argv[1]
    out_reduced = sys.argv[2]  # /home/anhvu/Desktop/mykointrons-data/reduced/
    scaffolds_portion = float(sys.argv[3])  # 0.04

    with open(fungi_names_file, 'r') as f:
        fungi_names = f.readlines()
        fungi_names = map(lambda name: name.rstrip(), fungi_names)

    for fungi_name in fungi_names:
        try:
            extract_fungi_portion(fungi_name, scaffolds_portion, out_reduced)
        except FileNotFoundError:
            print(f'Fungi {fungi_name} not found. Skipping')


def extract_fungi_portion(
        fungi_name: str,
        scaffolds_portion: float,
        out_reduced: str
):
    assembly = f'{ASSEMBLIES_LOC}/{fungi_name}_AssemblyScaffolds.fasta'
    introns_fasta = f'{NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-introns.fasta'

    out_reduced_assembly = f'{out_reduced}/Assembly/{fungi_name}_AssemblyScaffolds.fasta'
    out_reduced_introns = f'{out_reduced}/new-sequences/{fungi_name}-introns.fasta'

    with open(assembly, 'r') as assembly_f:
        # Load scaffolds in assembly
        scaffolds = SeqIO.parse(assembly_f, 'fasta')

        # Sample a set of assemblies
        scaffolds = sorted(scaffolds, key=lambda scaff: len(scaff.seq), reverse=True)
        scaffs_cumsum_lens = np.cumsum(list(map(lambda scaff: len(scaff.seq), scaffolds)))
        scaff_total_lens = scaffs_cumsum_lens[-1]

        scaffs_take_idx = np.argwhere(scaffs_cumsum_lens >= scaff_total_lens * scaffolds_portion)[0][0]

        print(f'{scaffs_take_idx + 1}/{len(scaffolds)} scaffolds will be taken from {fungi_name} '
              f'({100 * scaffs_cumsum_lens[scaffs_take_idx] / scaff_total_lens:.2f}% of genome)')

        taken_scaffolds = scaffolds[0:scaffs_take_idx + 1]
        taken_names = list(map(lambda scaff: scaff.id, taken_scaffolds))

        with open(out_reduced_assembly, 'w') as f:
            SeqIO.write(taken_scaffolds, f, 'fasta')

    with open(introns_fasta, 'r') as introns_f:
        introns = SeqIO.parse(introns_f, 'fasta')
        taken_introns = list(filter(lambda intron: intron.id in taken_names, introns))

        with open(out_reduced_introns, 'w') as f:
            SeqIO.write(taken_introns, f, 'fasta')


if __name__ == '__main__':
    main2()
