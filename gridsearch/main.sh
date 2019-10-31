#!/bin/bash

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
# Shuffle the examples and split them to train and validation sets.
# The train set will be used for grid search training

shroom_csv=../data/train/donor/Armme1_1-donor-windows.csv
bash split_csv_valid_train.sh $shroom_csv
# this will produce files with prefix 'shuffeled_train_' and 'shuffeled_valid_'

# Now perform grid search (on metacentre)
site=donor

bash ./grid_search $site "shuffeled_train_$shroom_csv"
# this will produce model files. Download them
model_folder=${site}_gridmodels
mkdir -p ./$model_folder

scp lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh/results/train-splice-site-shuffeled_train_Armme1_1-${site}* $model_folder

bash validate_grid_models.sh $site "shuffeled_valid_$shroom_csv" $model_folder

