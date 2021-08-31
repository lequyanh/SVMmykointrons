#!/bin/bash

# SYNOPSIS
#     bash combine_results.sh -p project_path
#
# OPTIONS
#     -p    project path with the fasta to be sharded into fragments (stored in ./assembly_shards)
# EXAMPLES
#     bash combine_results.sh -p /home/johndoe/Desktop/project/

while getopts "p:" opt; do
  case $opt in
  p)
    project_path=$OPTARG
    echo "Assembly FASTA for sharding will be taken from location ${project_path}"
    ;;
  *)
    echo "Invalid option or argument"
    exit 1
    ;;
  esac
done

results_folder="$project_path/results"

if [ -z "$results_folder" ]; then
  echo 'Result directory not found or empty. Exiting'
  exit
fi

assembly_path=$(find "$project_path" -maxdepth 1 -name "*no_duplicates*.fasta" -o -name "*.fa")
assembly_name=$(basename "$assembly_path" | cut -d '.' -f 1)

cut_coords_full_plus="$project_path/cut-coords-plus-strand-full.csv"
cut_coords_full_minus="$project_path/cut-coords-minus-strand-full.csv"

pruned_assembly_plus="$project_path/pruned_${assembly_name}_plus.fa"
pruned_assembly_minus="$project_path/pruned_${assembly_name}_minus.fa"

# ASSEMBLE PLUS STRAND RESULTS IF PRESENT
match="$(find "$project_path" -path "*results_plus*" -printf "." | wc -c)"
if [ "${match}" -gt 0 ]
then
  echo "scaffold;start;end" >"${cut_coords_full_plus}"

  for f in "${results_folder}"/*_results_plus/; do
    # Read all lines except for the first one containing headers
    tail -n +2 "$f/cut-coords.csv" >>"${cut_coords_full_plus}"
  done

  cat "${results_folder}"/*_results_plus/pruned-*.fa* > "${pruned_assembly_plus}"

  echo "Cleaned positive strand assembly completed"
fi

# ASSEMBLE MINUS STRAND RESULTS IF PRESENT
match="$(find "$project_path" -path "*results_minus*" -printf "." | wc -c)"

if [ "${match}" -gt 0 ]
then
  echo "scaffold;start;end" >"${cut_coords_full_minus}"

  for f in "${results_folder}"/*_results_minus/; do
    # Read all lines except for the first one containing headers
    tail -n +2 "$f/cut-coords.csv" >>"${cut_coords_full_minus}"
  done

  cat "${results_folder}"/*_results_minus/pruned-*.fa* > "${pruned_assembly_minus}"

  echo "Cleaned negative strand assembly completed"
fi
