#!/usr/bin/env bash

validation_set=$1
model_folders=$2

# Comment if extraction is not needed
for z in "${model_folders}"/*.tar.gz; do tar -xf "$z" --directory "$model_folders"; done
rm "${model_folders}"/*.tar.gz

for f in "${model_folders}"/*/
do
    order=$(basename "${f}" | cut -d "-" -f 6) 
    echo "Validation of model ${f} of order ${order}"
    qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=4gb:scratch_local=2gb -v MODEL_FOLDER="${f}",DATAFILE="${validation_set}",ORDER="${order}",CPU=16 validate_gridmodel_intron_meta.sh
done

