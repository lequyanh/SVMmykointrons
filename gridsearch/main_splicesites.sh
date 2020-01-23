#!/bin/bash

# ----- CHOOSE A SHROOM FOR GRID SEARCH TO INFER THE BEST SET OF SVM PARAMETERS FOR ACTUAL TRAINING ------
ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"

SITE="acceptor"
DIVISION="ascomycota"
SPECIES="Blade1"
SPECIES_CSV="../data/${DIVISION}/train/${SITE}/${SPECIES}-${SITE}-windows.csv"

TRAINSET="shuffeled_train_$(basename ${SPECIES_CSV})"
VALIDATION_SET="shuffeled_valid_$(basename ${SPECIES_CSV})"

# ratio of positive splice sites vs negative (false) splice_sites in the fungi genome
C_pos_real=0.003
# How to determine this parameter, go to:
#     ../statistics/splice-SITE-performance-evaluation/adjusting-precision-methodology.md
# But in general, for Ascomycota its roughly 0.003 and for Basidiomycota 0.025

# ======================================= GRID SEARCH FOR BEST SVM PARAMETERS =============================
function split_csv_valid_train() {
  echo "Shuffle the examples of chosen fungi species and split them to train and validation sets."
  echo "The train set will be used for grid search training"
  echo "Make sure there are enough examples (200k+)"

  bash split_csv_valid_train.sh ${SPECIES_CSV} "${TRAINSET}" "${VALIDATION_SET}"
}

function upload_trainset_and_train_gridmodels() {
  echo "Uploading train set to metacenter..."
  # Here we upload only the train file as training is more demanding and is better to do it remotely
  scp "${TRAINSET}" "${ROOT}/data/${SITE}s"

  echo "Now perform grid search (on metacenter). Use following command:"
  echo ">bash grid_search.sh ${SITE} ${TRAINSET}"
  echo "Change the range of parameters like C, win_in in the 'grid_search.sh' script according to your requirements"
}

function download_grid_models() {
  # local storage of grid models
  models_folder="${DIVISION}/${SITE}/models"
  mkdir -p "./${models_folder}"

  echo "Download trained grid models from metacentrum. Models will be save to ${models_folder}"

  # The result files contain the training filename in them
  scp "${ROOT}/results/${TRAINSET}*" ${models_folder}
}

function evaluate_grid_models() {
  results_folder="${DIVISION}/${SITE}/validation_results"
  mkdir -p "./${results_folder}"

  echo "Evaluate the performance of grid search models on the validation dataset."
  echo "Results will be stored in ${results_folder}"

  bash validate_grid_models.sh ${SITE} "${VALIDATION_SET}" ${models_folder} ${results_folder} ${C_pos_real}

  echo "Now we can compare the model performances to choose the best set of WD kernel parameters"
  echo "Guide on how to do so is in: ../statistics/splice-SITE-performance-evaluation/examining-gridsearch-result.md"

  # ------------------------------- LEGACY -----------------------------------
  # * The purpose of this is only to determine the adjusted precision. This metric is now already calculated during
  # validation of the model.
  # * Older result files however still need to be processed post-hoc, therefore this code stays

  # ratio of #true_splice_site_windows vs #false_splice_site_windows in resampled train/test file
  pos_ex=$(grep -c ";+1" <"${VALIDATION_SET}")
  neg_ex=$(grep -c ";-1" <"${VALIDATION_SET}")
  C_pos_resampled=$(awk "BEGIN {print $pos_ex/$neg_ex}")

  bash ../classification/evaluate_models.sh -f "${results_folder}" -r ${C_pos_real} -t "${C_pos_resampled}" -a 0.5 -g
  # ------------------------------- LEGACY -----------------------------------
}

split_csv_valid_train
upload_trainset_and_train_gridmodels
# Wait several hours...
#download_grid_models
#evaluate_grid_models

# =========================== TRAINING AND VALIDATION OF THE FULL MODEL WITH DISCOVERED PARAMETERS =====================
win_in=70
win_out=70
degree=25
C=10
model_name="${SITE}-model.hd5"

function train_full_model() {
  train_file="../data/${DIVISION}/train/shuff_aggreg_${SITE}_site_train.csv"

  if [ "$1" == "locally" ]; then
    python ../classification/train-splice-sites.py ${train_file} "${model_name}" ${win_in} ${win_out} ${degree} ${C} 12
  else
    echo "Upload the full train set to metacentrum"
    scp ${train_file} "${ROOT}/data/${SITE}s"

    echo "Execute this command on metacentrum (care for the interpratation of windows):"
    echo "qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=${degree},LWINDOW=70,RWINDOW=70,C=${C},SITE=${SITE}s,DATAFILE=$(basename ${train_file}),CPU=10 train-splice-SITE.sh"
    ## retrieve the model
  fi
}

function validate_full_model() {
  test_file="../data/${DIVISION}/test/shuff_aggreg_${SITE}_site_test.csv"
  C_pos_real=$1

  python ../classification/classify-splice-sites.py ${test_file} ${model_name} ${win_in} ${win_out} ${SITE} -r "${C_pos_real}" -c 10
}

#train_full_model "locally"
## For validation consider changing the C_+ parameter.
## But we assume that the ratio is more or less similar accross species
#validate_full_model "${C_pos_real}"
