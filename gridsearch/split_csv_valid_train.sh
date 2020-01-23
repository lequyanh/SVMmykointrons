#!/bin/bash

CSV_FILE=$1
TRAIN_CSV_FILE=$2
VALIDATION_CSV_FILE=$3

lines=$(wc -l < "${CSV_FILE}")
half=$((lines / 2))

# 400k is a reasonable size for SVM training - more would take too long to train
if [ $half -gt 400000 ] ; then
  tr_examples=400000
else
  tr_examples=$half
fi

# Read all except for the header (i.e. sequence;label)
tail -n +2 "${CSV_FILE}" > temp_seqs.txt
shuf temp_seqs.txt -o temp_seqs.txt


echo "sequence;label" > "${TRAIN_CSV_FILE}"
echo "sequence;label" > "${VALIDATION_CSV_FILE}"

{
  head -n $tr_examples >> "${TRAIN_CSV_FILE}"
  cat >> "${VALIDATION_CSV_FILE}"
} < temp_seqs.txt

rm temp_seqs.txt
