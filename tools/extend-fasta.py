import fileinput as fi
import fastalib as fl

'''
Accepts the input either from stdin or from a file (given by its filename as argument).
Prints the output to stdout.
'''

WIDTH = 70


def print_scaffold(desc, seq):
    print(">", desc, sep='')
    for i in range(0, len(seq), WIDTH):
        print(seq[i:i + WIDTH])


for description, sequence in fl.read_fasta(fi.input()):
    print_scaffold(description, sequence)

    negative_description = "neg___" + description
    negative_sequence = fl.complementary(reversed(sequence))

    print_scaffold(negative_description, negative_sequence)
