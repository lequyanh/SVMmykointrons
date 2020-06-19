#!/bin/bash

fmodel=/storage/praha1/home/lequyanh/data/ascomycota/models
dataf=/storage/praha1/home/lequyanh/data/reduced

dmodel="${fmodel}/donor-model-C-10-d-25-70.hd5"
amodel="${fmodel}/acceptor-model-C-10-d-25-win-70-shuffle2.hd5"
imodel="${fmodel}/intron-model-C-25-d-6.hd5"
phylum="ascomycota"

while read fungi; do
  assembly="${dataf}/Assembly/${fungi}_AssemblyScaffolds.fasta"
  introns="${dataf}/new-sequences/${fungi}-introns.fasta"
  echo "${fungi} with assembly ${assembly} and introns ${introns}"

  qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model=${dmodel},splice_site_acceptor_model=${amodel},intron_model=${imodel},phylum=${phylum},intron_source="${introns}",fungi_name="${fungi}" pipeline_perun.sh
done < "ascomycota.txt"
