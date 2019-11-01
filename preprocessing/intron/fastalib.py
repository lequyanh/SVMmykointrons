def read_fasta(file):
    '''
    Reads the input file and parses one scaffold at the time (memory efficient).
    '''
    description = None
    for line in file:
        if line.startswith('>'):
            if description is not None:
                yield description, sequence
            description = line.strip()[1:]
            sequence = ""
        else:
            sequence += line.strip()

    if sequence != "":
        yield description, sequence


complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def complementary(seq):
    '''
    Replaces nucleotides with their complements.
    '''
    return "".join(complements[c] if c in complements else c for c in seq)


def dimers(sequence):
    '''
    Generator for positions and dimers from the given sequence.
    '''
    return kmers(sequence, 2)


def kmers(sequence, k):
    '''
    Generator for positions and polymers from the given sequence.
    '''
    for i in range(len(sequence) - k + 1):
        yield (i + 1), sequence[i:i + k]
