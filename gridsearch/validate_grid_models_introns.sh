#!/usr/bin/env bash

NCPUS=12

validation_set=$1
model_folders=$2

# Comment if extraction is not needed
for z in "${model_folders}"/*.tar.gz; do tar -xf "$z" --directory "$model_folders"; done
rm "${model_folders}"/*.tar.gz

for f in "${model_folders}"/*/
do
    model="${f}/model.hd5"
    order=$(echo "$f" | cut -d "-" -f 6 )
    C_pos=$(echo "$f" | cut -d "-" -f 4 )

    echo "Evaluate intron model with d=${order} and C positive = ${C_pos}"

    python ../classification/classify-introns.py "${validation_set}" "${model}" "${order}" -c ${NCPUS}
    cat classify-introns*.log > "${f}/validation-results.txt"

    rm classify-introns*.log
done

