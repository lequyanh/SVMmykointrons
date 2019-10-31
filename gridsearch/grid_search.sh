#!/usr/bin/env bash

subject=$1
train_set=$2

ds=(15 20 25 30)
Cs=(1 10 100)
intron_windows=(60 70 80)

for deg in "${ds[@]}"; do
    for C in "${Cs[@]}"; do
        for i_win in "${intron_windows[@]}"; do
            if [ "$subject" == 'donor' ]
            then
                qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=${deg},LWINDOW=70,RWINDOW=${i_win},C=${C},SITE=donors,DATAFILE=${train_set},CPU=10 train-splice-site.sh
            else
                qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=${deg},LWINDOW=${i_win},RWINDOW=70,C=${C},SITE=acceptors,DATAFILE=${train_set},CPU=10 train-splice-site.sh
            fi
        done
    done
done
