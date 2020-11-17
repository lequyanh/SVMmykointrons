#!/bin/bash

# Parameters are:
# bash #this_script.sh -m #models_settings -l #list_of_fastas_to_process -d #directory_of_fastas

# Examples (validation):
# bash #this_script.sh -m svmb -l basidiomycota.txt -d /storage/praha1/home/lequyanh/data/reduced/ -s +
# Examples (classification on metagenom):
# bash #this_script.sh -m svmb -l metagenom_shards.txt -d /storage/praha1/home/lequyanh/data/metagenom/ -s + -f

while getopts "m:l:d:s:f" opt; do
  case $opt in
  m)
    models_settings=$OPTARG
    echo "Settings used for classification models ${models_settings}"
    ;;
  l)
    fasta_list_file=$OPTARG
    echo "The list of assemblies will be taken from ${fasta_list_file}"
    ;;
  d)
    fasta_dir=$OPTARG
    echo "Taking FASTAs from ${fasta_dir}"
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

if [ "$models_settings" == 'svma' ]; then
  fmodel="/storage/praha1/home/lequyanh/data/ascomycota/models"

  dmodel="${fmodel}/donor-model-C-10-d-25-70.hd5"
  amodel="${fmodel}/acceptor-model-C-10-d-25-win-70-shuffle2.hd5"
  imodel="${fmodel}/intron-model-C-25-d-6.hd5"
  window_inner=70
  window_outer=70

elif [ "$models_settings" == 'svmb' ]; then
  fmodel="/storage/praha1/home/lequyanh/data/basidiomycota/models"

  dmodel="${fmodel}/basidio-dmodel-plusminus-C-10-d-25-70.hd5"
  amodel="${fmodel}/basidio-amodel-plusminus-C-10-d-25-70-subsampled.hd5"
  imodel="${fmodel}/intron-model-C-7-d-6.hd5"
  window_inner=70
  window_outer=70

elif [ "$models_settings" == 'nnf' ]; then
  fmodel="/storage/praha1/home/lequyanh/data/nn_models"
  fmodel_basi="/storage/praha1/home/lequyanh/data/basidiomycota/models"

  dmodel="${fmodel}/bmodel.h5"
  amodel="${fmodel}/amodel.h5"
  imodel="${fmodel_basi}/intron-model-C-7-d-6.hd5"
  window_inner=200
  window_outer=200

elif [ "$models_settings" == 'nn100' ]; then
  fmodel="/storage/praha1/home/lequyanh/data/nn_models"
  fmodel_basi="/storage/praha1/home/lequyanh/data/basidiomycota/models"

  dmodel="${fmodel}/model_donor_100.h5"
  amodel="${fmodel}/model_acceptor_100.h5"
  imodel="${fmodel_basi}/intron-model-C-7-d-6.hd5"
  window_inner=100
  window_outer=0

elif [ "$models_settings" == 'random' ]; then
  dmodel="random"
  amodel="random"
  imodel="random"
  window_inner=70
  window_outer=70
fi

if [ $is_metagenom ]; then
  while read shard; do
    assembly="${fasta_dir}/${shard}"
    echo "Processing metagenom shard ${shard}"

    qsub -l walltime=48:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model="${dmodel}",splice_site_acceptor_model="${amodel}",intron_model="${imodel}",window_inner="${window_inner}",window_outer="${window_outer}",strand="${strand}",intron_source=None,fungi_name="${shard}" pipeline_perun.sh
  done <"$fasta_list_file"

else
  while read fungi; do
    assembly="${fasta_dir}/Assembly/${fungi}_AssemblyScaffolds.fasta"
    introns="${fasta_dir}/new-sequences/${fungi}-introns.fasta"
    echo "${fungi} with assembly ${assembly} and introns ${introns}"

    qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model="${dmodel}",splice_site_acceptor_model="${amodel}",intron_model="${imodel}",window_inner="${window_inner}",window_outer="${window_outer}",strand="${strand}",intron_source="${introns}",fungi_name="${fungi}" pipeline_perun.sh
  done <"$fasta_list_file"
fi
