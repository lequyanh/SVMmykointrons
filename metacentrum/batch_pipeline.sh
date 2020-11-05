#!/bin/bash

# Parameters are:
# bash #this_script.sh -p #phylum -l #list_of_fastas_to_process -d #directory_of_fastas

# Examples (validation):
# bash #this_script.sh -p basidiomycota -l basidiomycota.txt -d /storage/praha1/home/lequyanh/data/reduced/
# Examples (classification on metagenom):
# bash #this_script.sh -p basidiomycota -l metagenom_shards.txt -d /storage/praha1/home/lequyanh/data/metagenom/ -m


while getopts "p:l:d:m" opt; do
  case $opt in
  p)
    phylum=$OPTARG
    echo "The FASTAs will be processed by ${phylum} models"
    ;;
  l)
    fasta_list_file=$OPTARG
    echo "The list of assemblies will be taken from ${fasta_list_file}"
    ;;
  d)
    fasta_dir=$OPTARG
    echo "Taking FASTAs from ${fasta_dir}"
    ;;
  m)
    is_metagenom=TRUE
    echo "Taking FASTAs from ${fasta_dir}"
    ;;
  *)
    echo "Invalid option or argument"
    exit 1
    ;;
  esac
done

fmodel="/storage/praha1/home/lequyanh/data/${phylum}/models"

if [ "$phylum" == 'ascomycota' ]; then
  dmodel="${fmodel}/donor-model-C-10-d-25-70.hd5"
  amodel="${fmodel}/acceptor-model-C-10-d-25-win-70-shuffle2.hd5"
  imodel="${fmodel}/intron-model-C-25-d-6.hd5"

elif [ "$phylum" == 'basidiomycota' ]; then
  dmodel="${fmodel}/basidio-dmodel-plusminus-C-10-d-25-70.hd5"
  amodel="${fmodel}/basidio-amodel-plusminus-C-10-d-25-70-subsampled.hd5"
  imodel="${fmodel}/intron-model-C-7-d-6.hd5"

elif [ "$phylum" == 'random' ]; then
  dmodel="random"
  amodel="random"
  imodel="random"
fi

if [ $is_metagenom ]; then
  while read shard; do
    assembly="${fasta_dir}/${shard}"
    echo "Processing metagenom shard ${shard}"

    qsub -l walltime=48:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model="${dmodel}",splice_site_acceptor_model="${amodel}",intron_model="${imodel}",phylum="${phylum}",intron_source=None,fungi_name="${shard}" pipeline_perun.sh
  done <"$fasta_list_file"

else
  while read fungi; do
    assembly="${fasta_dir}/Assembly/${fungi}_AssemblyScaffolds.fasta"
    introns="${fasta_dir}/new-sequences/${fungi}-introns.fasta"
    echo "${fungi} with assembly ${assembly} and introns ${introns}"

    qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model="${dmodel}",splice_site_acceptor_model="${amodel}",intron_model="${imodel}",phylum="${phylum}",intron_source="${introns}",fungi_name="${fungi}" pipeline_perun.sh
  done <"$fasta_list_file"
fi
