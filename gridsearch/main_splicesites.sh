#!/bin/bash

ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
# Shuffle the examples and split them to train and validation sets.
# The train set will be used for grid search training
site="acceptor"
division="ascomycota"
shroom="Blade1"
shroom_csv=../data/${division}/train/${site}/${shroom}-${site}-windows.csv

# ======================================= GRID SEARCH FOR BEST SVM PARAMETERS =============================
bash split_csv_valid_train.sh $shroom_csv

# The command produced files with prefix 'shuffeled_train_' and 'shuffeled_valid_'. Upload them to metacenter
# Here we upload only the train file as training is more demanding and is better to do it remotely
scp shuffeled_train_${shroom}-${site}-windows.csv ${ROOT}/data/${site}s

echo "Now perform grid search (on metacenter). Use following command:"
echo "bash grid_search.sh $site shuffeled_train_$(basename ${shroom_csv})"

# ---------------------------------------------- Wait... --------------------------------------------------

# this will produce model files, which will be downloaded
models_folder=${site}_gridmodels/${division}
mkdir -p ./$models_folder
scp ${ROOT}/results/train-splice-site-shuffeled_train_${shroom}-${site}* $models_folder

# and evaluate their performance on validation dataset. Note the parameters of the best performing one
bash validate_grid_models.sh $site "shuffeled_valid_$(basename ${shroom_csv})" $models_folder "grid_validation_results/${division}"

# Compute metrics such as recall and adjusted precision. Change the parameters according to your data
C_pos_real=0.0033 # actual ratio of #true_splice_site_windows vs #false_splice_site_windows
# ratio of #true_splice_site_windows vs #false_splice_site_windows in resampled train/test file
# can be calculated from: cat shuffeled_train_${shroom}-{site}-windows.csv | grep ";-1" | wc -l
C_pos_resampled=0.25
bash ../classification/evaluate_models.sh -f "grid_validation_results/${division}/${site}/" -r ${C_pos_real} -t ${C_pos_resampled} -a 0.5 -g
# Examine the ouput to discover the best parameters

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
#echo "2) REMOTE: qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=30,LWINDOW=70,RWINDOW=70,C=1,SITE=donors,DATAFILE=$(basename ${train_file}),CPU=10 train-splice-site.sh"
## retrieve the model
## ------------------- metacenter ------------------

# ================================= VALIDATION OF THE MODEL =========================
# validate the model using the mixed shuffeled dataset (created during preprocessing step)
test_file="../data/${division}/test/shuff_aggreg_${site}_site_test.csv"
C_pos_resampled_test=0.33 # Ration of positive:negative examples in the validation set
model_location="../pipeline/bestmodels/${division}/${site}-model.hd5"

python ../classification/classify-splice-sites.py ${test_file} ${model_location} ${win_in} ${win_out} ${site} -v -c 10

result_file="overal-${site}-model-results.txt"

cat "classify-splice-sites.log" > ${result_file}
echo "Accuracy metrics recalculated with C_pos_real ${C_pos_real} and C_pos_resampled ${C_pos_resampled_test}" >> ${result_file}
bash ../classification/evaluate_models.sh -f "classify-splice-sites.log" -r ${C_pos_real} -t ${C_pos_resampled_test} -a 0.5 >> ${result_file}

rm "classify-splice-sites.log"