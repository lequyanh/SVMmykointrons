#!/bin/bash

DATA_LOC=$1   # Folder with results
OUT=$2        # Name of the output CSV file

echo "scaffold;start;end" > ${OUT}

while read shard; do
  # Read sequences in file
  # Read all lines except for the first one containing headers
  tail -n +2 "${DATA_LOC}/${shard}_results/cut-coords.csv" >> ${OUT}
done < metagenom_shards.txt

