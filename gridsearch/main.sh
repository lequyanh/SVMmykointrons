#!/bin/bash

ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
# Shuffle the examples and split them to train and validation sets.
# The train set will be used for grid search training
site="acceptor"
division="ascomycota"
shroom="Blade1"
shroom_csv=../data/${division}/train/${site}/${shroom}-${site}-windows.csv

bash split_csv_valid_train.sh $shroom_csv
# this will produce files with prefix 'shuffeled_train_' and 'shuffeled_valid_'. Upload them to metacenter
scp shuffeled_train_${shroom}-${site}-windows.csv ${ROOT}/data/${site}s

echo "Now perform grid search (on metacenter). Use following command:"
echo "bash ./grid_search $site shuffeled_train_$(basename ${shroom_csv})"

# Wait...

# this will produce model files, which will be downloaded
models_folder=${site}_gridmodels/${division}
mkdir -p ./$models_folder
scp ${ROOT}/results/train-splice-site-shuffeled_train_${shroom}-${site}* $models_folder

# and evaluate their performance on validation dataset. Note the parameters of the best performing one
bash validate_grid_models.sh $site "shuffeled_valid_$(basename ${shroom_csv})" $models_folder "grid_validation_results/${division}"
# Compute metrics such as recall and adjusted precision. Change the parameters according to your data
bash ../classification/evaluate_models.sh -f "grid_validation_results/${division}/${site}/" -r 0.025 -t 0.25 -a 0.5 -g

# train the overall model with the discovered parameters (locally or on metacenter)
# ------------------- locally ------------------
train_file="../data/${division}/train/shuff_aggreg_${site}_site_train.csv"
# Use the best parameters revealed by the grid search.
python ../classification/train-splice-sites.py ${train_file} "${site}-model.hd5" 70 70 30 1 10
 ------------------- metacenter ------------------
scp ${train_file} ${ROOT}/data/${site}s
echo "Use the following command to train the model:"
echo "qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=30,LWINDOW=70,RWINDOW=70,C=1,SITE=donors,DATAFILE=$(basename ${train_file}),CPU=10 train-splice-site.sh"
# retrieve the model
 ------------------- metacenter ------------------

# validate the model using the mixed shuffeled dataset (created during preprocessing step)
test_file="../data/${division}/test/shuff_aggreg_${site}_site_test.csv"
python ../classification/classify-splice-sites.py ${test_file} "model.hd5" 70 70 ${site} -v -c 10
#scp ../data/test/${test_file} lequyanh@skirit.metacentrum.cz:/~
#qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v VALIDATION_CSV=shuff_aggreg_donor_site_test.csv,MODEL=model.hd5,WIN_IN=70,SUBJECT=$site classify-splice-site.sh
