#!/bin/bash

# SYNOPSIS
#     bash shard_assembly.sh -p project_path -n seqs_per_shard
#
# OPTIONS
#     -p    project path with the fasta to be sharded into fragments (stored in ./assembly_shards)
#     -n    specify number of sequences per fasta shard (if not specified, defaults to 100)
# EXAMPLES
#     bash shard_assembly.sh -p /home/johndoe/Desktop/project/ -n 400

SEQ_PER_SHARD_DEFAULT=100
ASSEMBLY_SHARDS_DIR='assembly_shards'

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
      echo "Using default number of sequences per assembly shard $SEQ_PER_SHARD_DEFAULT"
      seqs_per_file=$SEQ_PER_SHARD_DEFAULT
fi

echo "Assembly will be split into shards with $seqs_per_file sequences each."

# Create a subdirectory in the project folder for assembly shards
assembly_shards_path="${project_path}/${ASSEMBLY_SHARDS_DIR}"
mkdir -p "${assembly_shards_path}"

if [ "$(ls -A "$assembly_shards_path")" ]; then
    echo "Directory with shards exists and is not empty: $assembly_shards_path"
    echo "Do you wish to overwrite the files? [y/n]"
    read -r yes_no

    if [ "$yes_no" != 'y' ]; then
        echo 'Aborting sharding and exiting'
        exit
    fi
fi

##############################
# SHARDING THE MAIN ASSEMBLY #
##############################
assembly_path=$(find "$project_path" -maxdepth 1 -name "*.fasta" -o -name "*.fa")
assembly_name=$(basename "$assembly_path" | cut -d '.' -f 1)

awk -v seqs_per_file="$seqs_per_file" \
    -v assembly_name="$assembly_name" \
    -v shards_dir="$assembly_shards_path" \
'BEGIN {n_seq=0;} /^>/ {
      if(n_seq % seqs_per_file == 0)
        { file=sprintf("%s/%s_%d.fa", shards_dir, assembly_name, n_seq); }
      print >> file;
      n_seq++;
      next;
    }
    {print >> file;}' < "$assembly_path"

echo "Assembly shards generated to the directory $assembly_shards_path"
