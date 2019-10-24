#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=30,LWINDOW=60,RWINDOW=80,C=1,SITE=donors,DATAFILE=data.csv.tar.gz,CPU=10 train-splice-site.sh

ROOT="/storage/praha1/home/$LOGNAME/"
PYTHONDIR="/storage/praha1/home/$LOGNAME/miniconda3/envs/mykointron/bin"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp $ROOT/data/$SITE/$DATAFILE $SCRATCHDIR || exit 1
cp $ROOT/scripts/train-splice-site.py $SCRATCHDIR || exit 1

cd $SCRATCHDIR || exit 2

# extract data
# tar -xzvf $DATAFILE

# select only some columns from the CSV - needs less memory afterwards
# cut -d';' -f4,5 *.csv > input.csv
cat *.csv > input.csv

MODELFILENAME="model.hd5"

$PYTHONDIR/python train-splice-site.py input.csv $MODELFILENAME $LWINDOW $RWINDOW $DEGREE $C $CPU > result.txt

mkdir -p $ROOT/results

RESULTDIR="train-splice-site-$DATAFILE-$C-$DEGREE-$LWINDOW-$RWINDOW-$PBS_JOBID"
mkdir $RESULTDIR

mv result.txt $RESULTDIR/
mv model.hd5 $RESULTDIR/

tar -czvf $RESULTDIR.tar.gz $RESULTDIR

cp $RESULTDIR.tar.gz $ROOT/results/

rm -rf $SCRATCHDIR/*
