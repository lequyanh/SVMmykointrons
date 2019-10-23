#!/usr/bin/env bash

# usage bash merge_shuffle_datasets.sh ~/PycharmProjects/mycointrons/data/train/donor 800000 shuff_aggreg_splice_site_train_800k.csv ~/PycharmProjects/mycointrons/data

DATA_LOC=$1
NO_EXAMPLES=$2
SHUFF_FILENAME=$3
OUTPUT_FOLDER=$4

TEMP_FOLDER='temp'
mkdir -p $TEMP_FOLDER

ls $DATA_LOC |sort -R |tail -20 | while read file; do
  cp $DATA_LOC/$file $TEMP_FOLDER
done

mkdir -p $TEMP_FOLDER/final
tail -n +2 $TEMP_FOLDER/*.csv > $TEMP_FOLDER/final/merged.csv
cd $TEMP_FOLDER/final

shuf -o shuffled.csv merged.csv
echo "sequence;label" > ${SHUFF_FILENAME}
head -n ${NO_EXAMPLES} shuffled.csv >> ${SHUFF_FILENAME}

#tar -czvf ${SHUFF_FILENAME}.tar.gz  ${SHUFF_FILENAME}
#cp ${SHUFF_FILENAME}.tar.gz /storage/praha1/home/lequyanh/data/

cd ../../
cp $TEMP_FOLDER/final/$SHUFF_FILENAME $OUTPUT_FOLDER
rm -r $TEMP_FOLDER

