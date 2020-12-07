#!/bin/bash

# Parameters are:
# bash #this_script.sh -m #models_settings -l #list_of_fastas_to_process -d #directory_of_fastas

# Examples (cutting and intron validation on fungal species):
# bash batch_pipeline.sh -m svmb -l basidiomycota.txt -d /storage/praha1/home/lequyanh/data/reduced/ -s +

# Examples (cutting on metagenom - no validation):
# bash batch_pipeline.sh -m svmb -l metagenom_shards.txt -d /storage/praha1/home/lequyanh/data/metagenom/ -s + -f

while getopts "m:l:d:s:f" opt; do
  case $opt in
  m)
    models_settings=$OPTARG
    echo "Settings used for classification models ${models_settings}"
    ;;
  l)
    fasta_list_file=$OPTARG
    echo "Names of assemblies to process will be taken from ${fasta_list_file}"
    ;;
  d)
    fasta_dir=$OPTARG
    echo "FASTAs of the assemblies will be taken from location ${fasta_dir}"
    ;;
  s)
    strand=$OPTARG
    echo "Strand ${strand}"
    ;;
  f)
    is_metagenom=TRUE
    echo "Flag to indicate we are working with metagenom set to True"
    ;;
  *)
    echo "Invalid option or argument"
    exit 1
    ;;
  esac
done

# Location of models
FMODEL_NN="/home/anhvu/PycharmProjects/mycointrons/pipeline/bestmodels/nn/"
FMODEL_BASI="/home/anhvu/PycharmProjects/mycointrons/pipeline/bestmodels/basidiomycota/"


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

#############################################################
# CLEANUP FUNCTION (prepares working directory for next run)#
#############################################################
function tidy(){
  assembly_name=$1

  target_dir="${assembly_name}_results"
  mkdir -p "$target_dir"

  mv ./*-result.csv "$target_dir"
  mv cut-coords.csv "$target_dir"
  mv pipeline.output "$target_dir"
  # mv intron-dataset.csv "$target_dir"

  rm ./*-dataset.csv
  rm pruned-*
}

###############################################
# META-GENOM CUTTING (FROM SHARDED FRAGMENTS) #
###############################################
if [ $is_metagenom ]; then
  while read shard; do
    assembly="${fasta_dir}/${shard}"
    echo "Processing metagenom shard ${shard}"

    bash pipeline "$assembly" $dmodel $amodel $imodel $window_inner $window_outer "$strand"
    tidy "$shard"
  done <"$fasta_list_file"

#################################
# FUNGI ASSEMBLY CUT VALIDATION #
#################################
else
  while read fungi; do
    assembly="${fasta_dir}/Assembly/${fungi}_AssemblyScaffolds.fasta"
    introns="${fasta_dir}/new-sequences/${fungi}-introns.fasta"
    echo "${fungi} with assembly ${assembly} and introns ${introns}"

    bash pipeline "$assembly" $dmodel $amodel $imodel $window_inner $window_outer "$strand" "$introns"
    tidy "$fungi"

  done <"$fasta_list_file"
fi
