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
                python ../classification/train-splice-sites.py "${train_set}" model.hd5 70 "${i_win}" "${deg}" "${C}" 12
                results_dir="train-${subject}-$(basename "$train_set")-C-${C}-d-${deg}-win-70-${i_win}"
            else
                python ../classification/train-splice-sites.py "${train_set}" model.hd5 "${i_win}" 70 "${deg}" "${C}" 12
                results_dir="train-${subject}-$(basename "$train_set")-C-${C}-d-${deg}-win-${i_win}-70"
            fi

            mkdir "${results_dir}"
            mv model.hd5 "${results_dir}"
            mv train-splice-sites-*.log "${results_dir}"
        done
    done
done
