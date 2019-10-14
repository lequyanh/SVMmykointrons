#!/usr/bin/env bash

# SCRATCHDIR='/c/Users/AnhVu/Desktop/scratch/'
# usage bash ggrid_search.sh ../data 500000 shrooms_all_shuffled_test.csv

DATA_LOC=$1
NO_EXAMPLES=$2
SHUFF_FILENAME=$3

cp ${DATA_LOC}/* $SCRATCHDIR
mkdir $SCRATCHDIR/final
cat $SCRATCHDIR/*.csv > $SCRATCHDIR/final/merged.csv
cd $SCRATCHDIR/final
shuf -o shuffled.csv merged.csv
head -n ${NO_EXAMPLES} shuffled.csv >  ${SHUFF_FILENAME}
tar -czvf ${SHUFF_FILENAME}.tar.gz  ${SHUFF_FILENAME}
cp ${SHUFF_FILENAME}.tar.gz /storage/praha1/home/lequyanh/data/
clean_scratch

