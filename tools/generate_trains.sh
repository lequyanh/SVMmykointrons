#!/bin/bash

PYTHON=/home/anhvu/anaconda3/envs/mykointron/bin/python

SHROOMS_FILE=$1
DATA_LOC=$2
WIN_IN=$3
WIN_OUT=$4
OUTPUT_FOLDER=$5
TRAIN_TEST=$6

for shroom in $(<$SHROOMS_FILE)
do
    $PYTHON ./extract-train-donor-acceptors.py $DATA_LOC $shroom $WIN_OUT $WIN_IN $WIN_IN $WIN_OUT $OUTPUT_FOLDER $TRAIN_TEST
done
