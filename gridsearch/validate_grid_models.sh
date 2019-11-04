#!/usr/bin/env bash

subject=$1
validation_set=$2
model_folders=$3
out_folder=$4

# Comment if extraction is not needed
for z in $model_folders/*.tar.gz; do tar -xf "$z" --directory $model_folders; done
rm $model_folders/*.tar.gz

result_dir="${out_folder}/${subject}"
mkdir -p "$result_dir"
echo "Grid validation in progress. Results will be saved in ./$result_dir/"

if [ "$subject" == "donor" ]; then
    col=10
else
    col=9
fi

for f in $model_folders/*/
do
    i_win=$(echo $f | cut -d'-' -f$col )
    model=$f/model.hd5

    python ../classification/classify-splice-sites.py $validation_set $model 70 $i_win $subject $result_dir -c 10
done