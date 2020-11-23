__author__ = 'vetrot'

import sys
import os
import operator

original_annotation = sys.argv[1]
fwd_annotation = sys.argv[2]
rev_annotation = sys.argv[3]
output_file = sys.argv[4]

class Annotation:
  def __init__(self, annotation):
      self.annotation = annotation
      self.bitscore = float(annotation.split('\t')[11])
      self.ctg = annotation.split('\t')[0].split('_')[1]
      self.A1 = int(annotation.split('\t')[0].split('_')[2])
      self.A2 = int(annotation.split('\t')[0].split('_')[3])
      self.intersections = []

  def getAnnotation(self):
    return self.annotation

  def getBestHit(self):
    self.intersections.sort(key=lambda x: x.bitscore, reverse=True)
    return self.intersections[0]

  def getInfo(self):
    print('From ' + str(self.A1) + ' to ' + str(self.A2) + ' on ' + self.ctg)

  def getIntersectInfo(self):
    print('Intersecting ' + str(len(self.intersections))+ ' anotations...')

  def getIntersectSize(self):
    return len(self.intersections)

  def addIntersectAnnot(self,b):
    self.intersections.append(b)

  def isIntersecting(self, b):
    index = 0
    intersect = False
    if self.A1 <= b.A1 and self.A2 >= b.A1:
        intersect = True
        index = 1
    # 2
    if b.A1 <= self.A1 and b.A2 >= self.A1:
        intersect = True
        index = 2
    # 3
    if self.A1 >= b.A1 and self.A2 <= b.A2:
        intersect = True
        index = 3
    # 4
    if self.A1 <= b.A1 and self.A2 >= b.A2:
        intersect = True
        index = 4
    #print(index)
    return intersect

def process_file(file, annot_rec):
    i = 0
    for line in open(file):
        a = Annotation(line.rstrip())
        if annot_rec.has_key(a.ctg):
            annot_rec[a.ctg].append(a)
        else:
            annot_rec[a.ctg] = []
            annot_rec[a.ctg].append(a)
        i += 1
    print("Anotaion " + file + " processed " + str(i) + " records...")


orig_annots = {}
process_file(original_annotation, orig_annots)

fwd_annots = {}
process_file(fwd_annotation, fwd_annots)

rev_annots = {}
process_file(rev_annotation, rev_annots)

for ctg in orig_annots:
    for a in orig_annots[ctg]:
        #a.getInfo()
        if fwd_annots.has_key(ctg):
            for b in fwd_annots[ctg]:
                #print('fwd '+str(a.isIntersecting(b))+'    [' + str(a.A1) + ',' + str(a.A2) + '] vs [' + str(b.A1)+','+str(b.A2)+']')
                if a.isIntersecting(b):
                    a.addIntersectAnnot(b)
                    b.addIntersectAnnot(a)
        if rev_annots.has_key(ctg):
            for b in rev_annots[ctg]:
                #print('rev '+str(a.isIntersecting(b))+'    [' + str(a.A1) + ',' + str(a.A2) + '] vs [' + str(b.A1)+','+str(b.A2)+']')
                if a.isIntersecting(b):
                    a.addIntersectAnnot(b)
                    b.addIntersectAnnot(a)
        #a.getIntersectInfo()

best_hits = open(output_file, "w")
replaced_hits = open(output_file+'_replaced.txt', "w")
new_hits = open(output_file+'_new.txt', "w")
# write best hit
io = 0
ir = 0
ip = 0
for ctg in orig_annots:
    for a in orig_annots[ctg]:
        io += 1
        if a.getIntersectSize() == 0:
            best_hits.write(a.getAnnotation() + '\n')
        else:
            b = a.getBestHit()
            if b.bitscore>a.bitscore:
                best_hits.write(b.getAnnotation() + '\n')
                replaced_hits.write(a.getAnnotation() + '\t' + b.getAnnotation() + '\n')
                ir += 1
            else:
                best_hits.write(a.getAnnotation() + '\n')
            # write possible best replacement...
            #replaced_hits.write(a.getAnnotation() + '\t' + b.getAnnotation() + '\n')
            ip += 1
print('Original hits ' + str(io) + ' replaced ' + str(ir)+ ' possible replacements ' + str(ip))

# compare fwd and reverse...
io = 0
ir = 0
ip = 0
for ctg in fwd_annots:
    for a in fwd_annots[ctg]:
        if a.getIntersectSize() == 0:
            if rev_annots.has_key(ctg):
                for b in rev_annots[ctg]:
                    if b.getIntersectSize() == 0:
                        if a.isIntersecting(b):
                            a.addIntersectAnnot(b)
                            b.addIntersectAnnot(a)
            # process it...
            if a.getIntersectSize() == 0:
                best_hits.write(a.getAnnotation() + '\n')
                new_hits.write(a.getAnnotation() + '\n')
                io += 1
            else:
                c = a.getBestHit()
                if c.bitscore>a.bitscore:
                    best_hits.write(c.getAnnotation() + '\n')
                    new_hits.write(c.getAnnotation() + '\n')
                    ir += 1
                else:
                    best_hits.write(a.getAnnotation() + '\n')
                    new_hits.write(a.getAnnotation() + '\n')
                ip += 1
print('FWD - not intersecting ' + str(io) + ' replaced ' + str(ir)+ ' possible replacements ' + str(ip))

# write not intersecting rev...
io = 0
for ctg in rev_annots:
    for b in rev_annots[ctg]:
        if b.getIntersectSize() == 0:
            best_hits.write(b.getAnnotation() + '\n')
            new_hits.write(b.getAnnotation() + '\n')
            io += 1
print('not intersecting rev ' + str(io))

best_hits.close()
replaced_hits.close()
new_hits.close()

print("DONE :]")

