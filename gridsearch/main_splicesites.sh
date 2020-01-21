#!/bin/bash

ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
# Shuffle the examples and split them to train and validation sets.
# The train set will be used for grid search training
site="acceptor"
division="ascomycota"
shroom="Blade1"
shroom_csv="../data/${division}/train/${site}/${shroom}-${site}-windows.csv"

# ======================================= GRID SEARCH FOR BEST SVM PARAMETERS =============================
bash split_csv_valid_train.sh $shroom_csv

# The command produced files with prefix 'shuffeled_train_' and 'shuffeled_valid_'. Upload them to metacenter
# Here we upload only the train file as training is more demanding and is better to do it remotely
scp "shuffeled_train_${shroom}-${site}-windows.csv" "${ROOT}/data/${site}s"

echo "Now perform grid search (on metacenter). Use following command:"
echo "bash grid_search.sh $site shuffeled_train_$(basename ${shroom_csv})"
echo "Change the range of parameters in the 'grid_search.sh' script according to your requirements"

# ---------------------------------------------- Wait... --------------------------------------------------

# this will produce model files, which will be downloaded
models_folder="${division}/${site}/models"
mkdir -p ./$models_folder
scp "${ROOT}/results/train-splice-site-shuffeled_train_${shroom}-${site}*" "$models_folder"

# and evaluate their performance on validation dataset. Note the parameters of the best performing one
shroom_validation_set="shuffeled_valid_$(basename ${shroom_csv})"
output_folder="${division}/${site}/validation_results"
mkdir -p ./$output_folder
bash validate_grid_models.sh $site "${shroom_validation_set}" "${models_folder}" "${output_folder}"

# ratio of #true_splice_site_windows vs #false_splice_site_windows in resampled train/test file
pos_ex=$(grep -c ";+1" < "${shroom_validation_set}")
neg_ex=$(grep -c ";-1" < "${shroom_validation_set}")
C_pos_resampled=$(awk "BEGIN {print $pos_ex/$neg_ex}")
# actual ratio of #true_splice_site_windows vs #false_splice_site_windows
C_pos_real=#CHANGE ME
# For details on how to set the two parameters, refer to:
# ../statistics/splice-site-performance-evaluation/adjusting-precision-methodology.md
bash ../classification/evaluate_models.sh -f "${output_folder}" -r ${C_pos_real} -t ${C_pos_resampled} -a 0.5 -g
# For more elaborate ways to choose hyperparameters from the grid-search results, please refer to:
# ../statistics/splice-site-performance-evaluation/examining-gridsearch-result.md

# ================================= TRAIN OVERALL MODEL WITH DISCOVERED PARAMETERS =========================
# Train the overall model with the discovered parameters (locally or on metacenter)
win_in=70
win_out=70
degree=30
C=1

# ------------------- locally ------------------
train_file="../data/${division}/train/shuff_aggreg_${site}_site_train.csv"
python ../classification/train-splice-sites.py ${train_file} "${site}-model.hd5" ${win_in} ${win_out} ${degree} ${C} 10
# ------------------- locally ------------------


## ------------------- metacenter ------------------
#echo "To train on metacenter, use the following commands:"
#echo "1) LOCALLY: scp ${train_file} ${ROOT}/data/${site}s"
#echo "2) REMOTLY: qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=30,LWINDOW=70,RWINDOW=70,C=1,SITE=donors,DATAFILE=$(basename ${train_file}),CPU=10 train-splice-site.sh"
## retrieve the model
## ------------------- metacenter ------------------

# ================================= VALIDATION OF THE MODEL =========================
# validate the model using the mixed shuffeled dataset (created during preprocessing step)
model_location="../pipeline/bestmodels/${division}/${site}-model.hd5"
test_file="../data/${division}/test/shuff_aggreg_${site}_site_test.csv"

pos_ex=$(grep -c ";+1" < "${test_file}")
neg_ex=$(grep -c ";-1" < "${test_file}")
C_pos_resampled_test=$(awk "BEGIN {print $pos_ex/$neg_ex}")

C_pos_real= # CHANGE_ME

python ../classification/classify-splice-sites.py ${test_file} ${model_location} ${win_in} ${win_out} ${site} -v -c 10

result_file="overal-${site}-model-results.txt"

cat "classify-splice-sites.log" > ${result_file}
echo "Accuracy metrics recalculated with C_pos_real ${C_pos_real} and C_pos_resampled ${C_pos_resampled_test}" >> ${result_file}
bash ../classification/evaluate_models.sh -f "classify-splice-sites.log" -r ${C_pos_real} -t "${C_pos_resampled_test}" -a 0.5 >> ${result_file}

rm "classify-splice-sites.log"