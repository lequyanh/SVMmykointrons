#!/bin/bash

aggreg_filename="$(pwd)/ascomycota-intron-lens.txt"
new_sequences_loc="/home/anhvu/Desktop/mykointrons-data/new-sequences/"

for f in ${new_sequences_loc}/*/; do
  cd "$f" || exit
  cat "$(basename ${f})-introns-len.txt" > $aggreg_filename
done