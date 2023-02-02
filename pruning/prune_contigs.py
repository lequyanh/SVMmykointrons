__author__ = 'vetrot'

import sys
import operator
import os

input_fasta = sys.argv[1]
coords_file = sys.argv[2]
pruned_fasta = sys.argv[3]

contigs = {}
reduction = {}
filled = False
for n, line1 in enumerate(open(input_fasta)):
    if n % 200000 == 0:
        print(str(n / 2))
    if n % 2 == 0:
        title = line1.rstrip()[1:].split(' ')[0]
    else:
        if n % 2 == 1:
            sequence = line1.rstrip()
            filled = True
    if filled:
        contigs[title] = sequence
        reduction[title] = 0
        filled = False

print("contigs read...")

for n, line in enumerate(open(coords_file)):
    if n == 0:
        header = line.rstrip()
    else:
        vals = line.rstrip().split(';')
        contig = contigs[vals[0]]
        start = int(vals[1]) - reduction[vals[0]]
        stop = int(vals[2]) - reduction[vals[0]]
        #print(vals[0] +' reduction: '+str(reduction[vals[0]])+' start: ' + str(start) +'('+vals[1]+')' + ' stop: ' + str(stop) +'(' + vals[2] + ')')
        reduction[vals[0]] += (stop - start) + 1
        # cut and store...
        contigs[vals[0]] = contig[:start - 1] + contig[stop:]

print("pruning done...")

fp = open(pruned_fasta, "w")
for title in contigs:
    #print(title)
    #print(contigs[title])
    fp.write(">" + title + "\n")
    fp.write(contigs[title] + "\n")

fp.close()

print("Output saved....")
