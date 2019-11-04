#!/bin/bash

ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
# Shuffle the examples and split them to train and validation sets.
# The train set will be used for grid search training
site="donor"

shroom="Armme1_1"
shroom_csv=../data/train/${site}/${shroom}-donor-windows.csv

bash split_csv_valid_train.sh $shroom_csv
# this will produce files with prefix 'shuffeled_train_' and 'shuffeled_valid_'. Upload them to metacenter
scp shuffeled_train_${shroom}-${site}-windows.csv ${ROOT}/data/${site}s

# Now perform grid search (on metacenter)
bash ./grid_search $site "shuffeled_train_${shroom_csv}"

# this will produce model files. Download them
models_folder=${site}_gridmodels
mkdir -p ./$models_folder
scp ${ROOT}/results/train-splice-site-shuffeled_train_${shroom}-${site}* $models_folder

# and evaluate their performance on validation dataset. Note the parameters of the best performing one
bash validate_grid_models.sh $site "shuffeled_valid_$shroom_csv" $models_folder

# train the overall model with the discovered parameters (on metacenter)
train_file="shuff_aggreg_${site}_site_train.csv"
scp ../data/train/${train_file} ${ROOT}/data/${site}s
qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=30,LWINDOW=70,RWINDOW=70,C=1,SITE=donors,DATAFILE=${train_file},CPU=10 train-splice-site.sh
# retrieve the model

# validate the model using the mixed shuffeled dataset (created during preprocessing step)
test_file="shuff_aggreg_${site}_site_test.csv"
scp ../data/test/${test_file} lequyanh@skirit.metacentrum.cz:/~
qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v VALIDATION_CSV=shuff_aggreg_donor_site_test.csv,MODEL=model.hd5,WIN_IN=70,SUBJECT=$site classify-splice-site.sh
