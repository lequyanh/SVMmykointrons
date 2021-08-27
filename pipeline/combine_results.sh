#!/bin/bash

# Function for merging cuts into one file
function merge_cuts() {
  results_folder=$1
  out=$2

  echo "scaffold;start;end" >"${out}"

  for f in "${results_folder}"/*/; do
    # Read all lines except for the first one containing headers
    tail -n +2 "$f/cut-coords.csv" >>"${out}"
  done
}
# Output file names
OUT_FORW='metagenom-cut-coords.csv'
OUT_REV='metagenom-cut-coords-reverse.csv'

merge_cuts ${forw_strand_result_loc} $OUT_FORW
merge_cuts ${rev_strand_result_loc} $OUT_REV
