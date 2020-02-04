#!/usr/bin/env bash

subject=$1
validation_set=$2
models_folders=$3

# Comment if extraction is not needed
for z in "${models_folders}"/*.tar.gz; do tar -xf "$z" --directory "${models_folders}"; done
rm "${models_folders}"/*.tar.gz

if [ "$subject" == "donor" ]; then
    col=10
else
    col=9
fi

for f in "${models_folders}"/*/
do
    i_win=$(echo "${f}" | cut -d'-' -f$col )

    qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v MODEL_FOLDER="$f",I_WIN="${i_win}",C=1,SUBJECT="${subject}",DATAFILE="${validation_set}",IMBALANCE_RAT=0.025,CPU=10 cloud_classify_gridmodel.sh

done

