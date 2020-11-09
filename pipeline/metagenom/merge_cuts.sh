#!/bin/bash

DATA_LOC=random_metagenom_results
OUT=metagenom-cut-coords-random.csv

echo "scaffold;start;end" > ${OUT}

while read shard; do
  # Read sequences in file
  # Read all lines except for the first one containing headers
  tail -n +2 "${DATA_LOC}/${shard}_results/cut-coords.csv" >> ${OUT}
done < metagenom_shards.txt

