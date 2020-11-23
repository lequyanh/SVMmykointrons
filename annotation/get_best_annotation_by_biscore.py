__author__ = 'vetrot'

import sys
import os

annotation1 = sys.argv[1]
annotation2 = sys.argv[2]
output_file = sys.argv[3]

annot = {}
for line in open(annotation1):
    vals = line.rstrip().split('\t')
    annot[vals[0]] = line.rstrip()

print('Annotation 1 was loaded - ' + str(len(annot)) + ' records')

i = 0
n = 0
fp = open(output_file, "w")
for line in open(annotation2):
    vals = line.rstrip().split('\t')
    if annot.has_key(vals[0]):
        #compare
        a1 = annot[vals[0]]
        bitscore = float(a1.split('\t')[11])
        #print("a1 score: "+str(bitscore)+" a2 score: "+vals[11]+" a1 eval: "+a1.split('\t')[10]+" a2 eval: "+vals[10])
        if bitscore > float(vals[11]):
            fp.write(a1 + '\n')
            n += 1
        else:
            fp.write(line.rstrip() + '\n')
        #clean
        annot[vals[0]] = ''
    else:
        fp.write(line.rstrip() + '\n')
    i += 1

print('Second annotation was compared and written - ' + str(i) + ' records - replaced ' + str(n))

for key in annot:
    if not annot[key] == '':
        fp.write(annot[key] + '\n')
        i += 1

print('Rest of annotations was written - total ' + str(i) + ' records')

fp.close()

print('Done :]')