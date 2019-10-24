#!/usr/bin/env bash

# usage bash merge_shuffle_datasets.sh ~/PycharmProjects/mycointrons/data/train/donor 800000 shuff_aggreg_splice_site_train_800k.csv ~/PycharmProjects/mycointrons/data

DATA_LOC=$1
NO_EXAMPLES=$2
SHUFF_FILENAME=$3
OUTPUT_FOLDER=$4

echo "Examples in [$DATA_LOC] will be merged and shuffled."
echo "[$NO_EXAMPLES] is the target number of examples in an aggregated train/test set"

TEMP_FOLDER='temp'
rm -r $TEMP_FOLDER
mkdir $TEMP_FOLDER

NO_REPR=$(( $NO_EXAMPLES/$(ls $DATA_LOC | wc -l) ))
echo "[$NO_REPR] random representatives will be chosen from each CSV file for merging a shuffling"

ls $DATA_LOC | while read file; do
  # Read sequences in file
  # Read all lines except for the first one containing headers
  tail -n +2 $DATA_LOC/$file > temp_seqs.txt
  shuf temp_seqs.txt | head -n $NO_REPR >> $TEMP_FOLDER/merged.csv
done

cd $TEMP_FOLDER/

shuf -o shuffled.csv merged.csv
echo "sequence;label" > ${SHUFF_FILENAME}
head -n ${NO_EXAMPLES} shuffled.csv >> ${SHUFF_FILENAME}

echo "Examples saved to [$SHUFF_FILENAME]"

#tar -czvf ${SHUFF_FILENAME}.tar.gz  ${SHUFF_FILENAME}
#cp ${SHUFF_FILENAME}.tar.gz /storage/praha1/home/lequyanh/data/

cd ../
cp $TEMP_FOLDER/$SHUFF_FILENAME $OUTPUT_FOLDER
rm -r $TEMP_FOLDER

