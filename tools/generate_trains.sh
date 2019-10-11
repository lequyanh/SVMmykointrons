#!/bin/bash

PYTHON=/home/anhvu/anaconda3/envs/mykointron/bin/python

SHROOMS_FILE=$1
WIN_IN=$2
WIN_OUT=$3

for shroom in $(<$SHROOMS_FILE)
do
   $PYTHON ./extract-train-donor-acceptors.py '/home/anhvu/Desktop/mykointrons-data' $shroom  $WIN_OUT $WIN_IN $WIN_IN $WIN_OUT '../data/splice-sites-train'
done
