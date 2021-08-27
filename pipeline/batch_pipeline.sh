#!/bin/bash

# Basic usage:
# bash batch_pipeline.sh -m #models_settings -p #project_path -s #strand [-l #custom_shards_list]

# EXAMPLES
# bash batch_pipeline.sh -m nn100 -p /home/johndoe/Desktop/project/ -s plus

# Specifying only a subset of assembly shards to process (using SVMs)
# bash main.sh -m svmb -p /home/johndoe/Desktop/project/ -l custom_assembly_shard_list.txt

# OBSOLETE: cutting and intron validation on fungal species):
# bash batch_pipeline.sh -m svmb -l basidiomycota.txt -d /storage/praha1/home/lequyanh/data/reduced/ -s minus

####################################################################
# FUNCTION FOR DISTRIBUTED COMPUTATION                             #
#    Please modify the function by inserting a job submission call #
#    according to standards of your cluster                        #
####################################################################
function submit_job(){
  assembly=$1
  working_dir=$2
  results_dir=$3

  # >>>> !! CHANGE ME !! <<<
  # Use this call for local execution or replace it a job submission call
  bash pipeline "$assembly" "$dmodel" "$amodel" "$imodel" "$window_inner" "$window_outer" "$strand" "$working_dir" "$results_dir"
  # >>>> !! CHANGE ME !! <<<
}
export -f submit_job

######################
# Defining constants #
######################
# Location of models and the scripts
ROOT=$(dirname "$(pwd)")
FMODEL_NN="$(pwd)/bestmodels/nn/"
FMODEL_BASI="$(pwd)/bestmodels/basidiomycota/"

# Project sub-folders for partial and full results
SHARDS_DIR="assembly_shards"
SCRATCH_DIR="scratchdir"
RESULTS_DIR="results"

export SHARDS_DIR
export SCRATCH_DIR
export RESULTS_DIR

########################
# Processing arguments #
########################
while getopts "m:p:s:l:" opt; do
  case $opt in
  m)
    models_settings=$OPTARG
    echo "Settings used for classification models ${models_settings}"
    ;;
  p)
    project_path=$OPTARG
    echo "FASTAs of the assemblies will be taken from location ${SHARDS_DIR}"
    ;;
  l)
    fasta_list_file=$OPTARG
    echo "Names of assemblies to process will be taken from ${fasta_list_file}"
    ;;
  s)
    strand=$OPTARG
    echo "Strand ${strand}"
    ;;
  *)
    echo "Invalid option or argument"
    exit 1
    ;;
  esac
done

if [ -z "$fasta_list_file" ]
then
      echo "Processing all assembly assembly_shards in $project_path"

      fasta_list_file="${project_path}/assembly_shards.txt"
      ls "${project_path}/${SHARDS_DIR}" >"${fasta_list_file}"
else
      echo "Using a subset of assembly assembly_shards defined in $fasta_list_file"
fi

export strand
export project_path

###############################################
# PIPELINE PARAMETER SETTINGS BASED ON MODEL  #
###############################################
if [ "$models_settings" == 'svmb' ]; then
  dmodel="${FMODEL_BASI}/basidio-dmodel-plusminus-C-10-d-25-70.hd5"
  amodel="${FMODEL_BASI}/basidio-amodel-plusminus-C-10-d-25-70-subsampled.hd5"
  imodel="${FMODEL_BASI}/intron-model-C-7-d-6.hd5"
  window_inner=70
  window_outer=70

elif [ "$models_settings" == 'nn200' ]; then
  dmodel="${FMODEL_NN}/model_donor_200.h5"
  amodel="${FMODEL_NN}/model_acceptor_200.h5"
  imodel="${FMODEL_BASI}/intron-model-C-7-d-6.hd5"
  window_inner=200
  window_outer=200

elif [ "$models_settings" == 'nn100' ]; then
  dmodel="${FMODEL_NN}/model_donor_100.h5"
  amodel="${FMODEL_NN}/model_acceptor_100.h5"
  imodel="${FMODEL_BASI}/intron-model-C-7-d-6.hd5"
  window_inner=100
  window_outer=0

elif [ "$models_settings" == 'random' ]; then
  dmodel="random"
  amodel="random"
  imodel="random"
  window_inner=70
  window_outer=70
fi

export amodel
export dmodel
export imodel
export window_inner
export window_outer

#############################################################
# FUNCTIONS (prepares working directory for next run)#
#############################################################
function process_assembly_assembly_shard(){
  assembly_shard=$1

  # Prepare working and results directories for independent processing
  assembly="${project_path}/${SHARDS_DIR}/${assembly_shard}"
  working_dir="${project_path}/${SCRATCH_DIR}/${assembly_shard}_${strand}"
  results_dir="${project_path}/${RESULTS_DIR}/${assembly_shard}_results_${strand}"

  mkdir -p "$working_dir"
  ln -s -f "$ROOT/classification/classify-splice-sites.py" "${working_dir}/classify-splice-sites.py"
  ln -s -f "$ROOT/classification/classify-introns.py" "${working_dir}/classify-introns.py"
  ln -s -f "$ROOT/classification/classify-introns-stub.py" "${working_dir}/classify-introns-stub.py"
  ln -s -f "$ROOT/classification/tools.py" "${working_dir}/tools.py"

  ln -s -f "$ROOT/tools/fastalib.py" "${working_dir}/fastalib.py"
  ln -s -f "$ROOT/tools/generate-pairs.py" "${working_dir}/generate-pairs.py"

  ln -s -f "$ROOT/pruning/prune_tools.py" "${working_dir}/prune_tools.py"
  ln -s -f "$ROOT/pruning/prune_probabilistic.py" "${working_dir}/prune_probabilistic.py"

  ln -s -f "$ROOT/pipeline/filter-orphan-acceptors.py" "${working_dir}/filter-orphan-acceptors.py"
  ln -s -f "$ROOT/pipeline/extract-donor-acceptor.py" "${working_dir}/extract-donor-acceptor.py"
  ln -s -f "$ROOT/pipeline/extract-introns.py" "${working_dir}/extract-introns.py"
  ln -s -f "$ROOT/pipeline/basidiomycota-intron-lens.txt" "${working_dir}"

  # tell a HPC to run the pipeline script with these parameters (register a job or something similar)
  submit_job "$assembly" "$working_dir" "$results_dir"
}

###############################################
# META-GENOM CUTTING (FROM SHARDED FRAGMENTS) #
###############################################
export -f process_assembly_assembly_shard

while read assembly_shard; do
  echo "Sending for execution the assembly shard ${assembly_shard}"
  process_assembly_assembly_shard "${assembly_shard}"
done <"$fasta_list_file"
