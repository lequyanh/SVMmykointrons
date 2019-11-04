#!/bin/bash

CSV_FILE=$1

lines=$(wc -l < "$CSV_FILE")
half=$((lines / 2))

# 400k is a reasonable size for SVM training - more would take too long to train
if [ $half -gt 400000 ] ; then
  tr_examples=400000
else
  tr_examples=$half
fi

# Read all except for the header (i.e. sequence;label)
tail -n +2 "$CSV_FILE" > temp_seqs.txt
shuf temp_seqs.txt -o temp_seqs.txt

train_file="shuffeled_train_$(basename "$CSV_FILE")"
validation_file="shuffeled_valid_$(basename "$CSV_FILE")"
echo "sequence;label" > "$train_file"
echo "sequence;label" > "$validation_file"

{
  head -n $tr_examples >> "$train_file"
  cat >> "$validation_file"
} < temp_seqs.txt

rm temp_seqs.txt
