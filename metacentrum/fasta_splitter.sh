#!/bin/bash

# Splits large FASTA (usually metagenom) into smaller shard, that can be then further processed in parallel

LARGE_FASTA=$1
NO_SHARDS=$2

lines=$(wc -l <"$LARGE_FASTA")
echo "Number of lines ${lines}. Will be split to cca ${NO_SHARDS} files".

lines_per_shard=$(("$lines" / "$NO_SHARDS"))
if [ $((lines_per_shard % 2)) == 1 ]; then
  lines_per_shard=$((lines_per_shard + 1))
fi

split -l "${lines_per_shard}" "${LARGE_FASTA}" --numeric-suffixes "metagenom_shard" --suffix-length=${#NO_SHARDS}
