#!/usr/bin/env bash

NCPUS=12

validation_set=$1
model_folders=$2
out_folder=$3

# Comment if extraction is not needed
for z in ${model_folders}/*.tar.gz; do tar -xf "$z" --directory "$model_folders"; done
rm "${model_folders}/*.tar.gz"

mkdir -p "$out_folder"
echo "Grid validation of models in ${model_folders}. Results will be saved in ./$out_folder"

for f in $model_folders/*/
do
    model=$f/model.hd5
    order=$(echo "$f" | cut -d "-" -f 4 )
    C_pos=$(echo "$f" | cut -d "-" -f 3 )

    echo "Evaluate intron model with d=${order} and C positive = ${C_pos}"

    # ---------- CLOUD VERSION --------------
    # qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v SUBJECT=$subject,DATA=$validation_set,MODEL=$model_folders,I_WIN=$i_win,CPU=10 classify-splice-sites-cloud.sh
    # ---------- CLOUD VERSION --------------

    # ---------- LOCAL VERSION --------------
    python ../classification/classify-introns.py "$validation_set" "$model" "$order" -c $NCPUS
    cat 'classify-introns.log' > "./${out_folder}/$(basename "$f")-results.txt"
    # ---------- LOCAL VERSION --------------
done

rm *.log