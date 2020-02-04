#!/usr/bin/env bash

subject=$1        # donor/acceptor
validation_set=$2
model_folders=$3  # typically "gridsearch/ascomycota/$subject/models"
result_dir=$4     # typically "gridsearch/ascomycota/$subject/validation_results"
pos_neg_ratio=$5  # ratio between positive and negative (false) splice sites in mushroom genome

# Comment if extraction is not needed
for z in "${model_folders}"/*.tar.gz; do tar -xf "$z" --directory "${model_folders}"; done
rm "${model_folders}"/*.tar.gz

mkdir -p "${result_dir}"
echo "Grid models validation in progress. Results will be saved in ./${result_dir}"

# Folder with models have parameters of that model in their name (split by '-').
# Extract the INNER_WINDOW parameter from the names
if [ "$subject" == "donor" ]; then
    col=12
else
    col=11
fi

for f in "${model_folders}"/*/
do
    i_win=$(echo "${f%?}" | cut -d'-' -f$col )  # the term after echo removes the last character
    model="${f}/model.hd5"

    python ../classification/classify-splice-sites.py "${validation_set}" "${model}" "${i_win}" 70 "${subject}" -r "${pos_neg_ratio}" -c 12
    cat classify-splice-sites*.log > "./${result_dir}/$(basename "$f")-results.txt"
done
