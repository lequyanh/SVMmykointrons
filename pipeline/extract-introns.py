import sys
from collections import defaultdict
from fastalib import read_fasta


def sequences(file):
    for line in file.readlines():
        words = line.strip().split(' ')
        assert len(words) >= 3,\
            'each line must contain at least scaffold with a start and stop position'
        assert len(words) % 2 == 1,\
            'each line must contain a scaffold and pairs of start/end positions'
        yield words[0], [int(p) for p in words[1:]]


if __name__ == '__main__':
    #  name of the FASTA assembly file
    fasta_filename = sys.argv[1]
    # name of the file containing the sequences
    sequences_filename = sys.argv[2]

    # scaffold start end start end start end

    scaffold_sequences = defaultdict(list)

    # load all sequences (a sequence is identified by the scaffold and start/end position pairs)
    with open(sequences_filename, 'r') as sequences_file:
        for scaffold, positions in sequences(sequences_file):
            scaffold_sequences[scaffold].append(positions)

    print('scaffold;start;end;sequence')
    with open(fasta_filename, 'r') as fasta_file:
        for scaffold, sequence in read_fasta(fasta_file):
            to_extract = scaffold_sequences[scaffold]

            for positions in to_extract:
                extracted = ''.join(sequence[positions[i] - 1:positions[i+1]]
                                    for i in range(0, len(positions), 2))
                print(scaffold, *positions, extracted, sep=';')

