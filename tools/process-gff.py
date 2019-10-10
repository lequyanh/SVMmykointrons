import fileinput as fi
from collections import defaultdict
from re import compile

'''
Accepts the input either from stdin or from a GFF file (given by its filename as argument).
Prints the output to stdout.
'''


def extract_introns(sequences):
    """
    Returns positions of introns in a gene given non-intron sequences
    present in the gene.
    """
    sequences = map(lambda s: (s[2], s[3]) if s[2] <= s[3] else (s[3], s[2]),
                    sequences)
    sequences = sorted(sequences, key=lambda s: s[0])

    s = list(map(lambda s: s[0], sequences))
    e = list(map(lambda s: s[1], sequences))

    start = min(s)
    end = max(e)

    last = start
    introns = []

    for i in range(len(sequences)):
        if s[i] > last + 1:
            introns.append((last + 1, s[i] - 1))
        last = max(last, e[i])

    return introns


name_regex = compile("^name \"(.*?)\"")

genes = defaultdict(list)

for line in fi.input():
    line = line.strip()
    cols = line.split('\t')

    scaffold, start, end, strand, attrs = \
        cols[0], cols[3], cols[4], cols[6], cols[8]

    start, end = int(start), int(end)
    gene_name = name_regex.match(attrs).group(1)

    genes[gene_name].append((scaffold, strand, start, end))

for gene, sequences in genes.items():
    introns = extract_introns(sequences)

    scaffold, strand = sequences[0][:2]

    for intron in introns:
        print(scaffold, strand, intron[0], intron[1])
