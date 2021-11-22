#!/bin/bash

# SYNOPSIS
#     bash shard_assembly.sh -p project_path -n seqs_per_batch
#
# OPTIONS
#     -p    project path with the fasta to be sharded into fragments (stored in ./assembly_batches)
#     -n    specify number of sequences per fasta batch (if not specified, defaults to 100)
# EXAMPLES
#     bash shard_assembly.sh -p /home/johndoe/Desktop/project/ -n 400

SEQ_PER_BATCH_DEFAULT=100
ASSEMBLY_BATCHES_DIR='assembly_batches'

while getopts "p:n:" opt; do
  case $opt in
  p)
    project_path=$OPTARG
    echo "Assembly FASTA for sharding will be taken from location ${project_path}"
    ;;
  n)
    seqs_per_file=$OPTARG
    ;;
  *)
    echo "Invalid option or argument"
    exit 1
    ;;
  esac
done

if [ -z "$seqs_per_file" ]
then
      echo "Using default number of sequences per assembly batch $SEQ_PER_BATCH_DEFAULT"
      seqs_per_file=$SEQ_PER_BATCH_DEFAULT
fi
echo "Assembly will be split into batches with $seqs_per_file sequences each."

no_fastas=$(find "$project_path" -maxdepth 1 -name "*.fasta" -o -name "*.fa" -not -name "*_no_duplicates.fa" | wc -l)
if [ "$no_fastas" -gt 1 ]; then
  echo "Multiple FASTA files in the project directory!"
  exit 1
fi

# Create a subdirectory in the project folder for assembly batches
assembly_batches_path="${project_path}/${ASSEMBLY_BATCHES_DIR}"
mkdir -p "${assembly_batches_path}"

if [ "$(ls -A "$assembly_batches_path")" ]; then
    echo "Directory with assembly batches exist and is not empty: $assembly_batches_path"
    echo "Do you wish to overwrite the files? [y/n]"
    read -r yes_no

    if [ "$yes_no" != 'y' ]; then
        echo 'Aborting sharding and exiting'
        exit
    fi

    rm "$assembly_batches_path"/*
fi

##############################
# SHARDING THE MAIN ASSEMBLY #
##############################
assembly_path=$(find "$project_path" -maxdepth 1 -name "*.fasta" -o -name "*.fa" -not -name "*_no_duplicates.fa")
assembly_name=$(basename "$assembly_path" | cut -d '.' -f 1)

# Remove duplicated entries (save to a separate fasta)
assembly_unique="${project_path}/${assembly_name}_no_duplicates.fa"

printf "\nRemoving duplicated sequences..."
if [ -z "${assembly_unique}" ]
then
  printf "\tNon-duplicated fasta existing. Proceed to sharding"
else
  printf "\tOriginal number of sequences: %s \n"  "$(grep -c '>' < "$assembly_path")"

  awk '/^>/ { f = !($0 in a); a[$0] } f' "$assembly_path" > "$assembly_unique"

  printf "\tNumber of unique sequences after cleaning: %s \n" "$(grep -c '>' < "$assembly_unique")"
  printf "Non-duplicated sequences are saved to %s \n\n" "$assembly_unique"
fi

# Sharding
echo "Sharding the assembly..."

awk -v seqs_per_file="$seqs_per_file" \
    -v assembly_name="$assembly_name" \
    -v batches_dir="$assembly_batches_path" \
'BEGIN {n_seq=0;} /^>/ {
      if(n_seq % seqs_per_file == 0)
        { file=sprintf("%s/%s_%d.fa", batches_dir, assembly_name, n_seq); }
      print >> file;
      n_seq++;
      next;
    }
    {print >> file;}' < "$assembly_unique"

echo "Assembly batches generated to the directory $assembly_batches_path"
