fmodel=/storage/praha1/home/lequyanh/data/basidiomycota/models
dataf=/storage/praha1/home/lequyanh/data/reduced

dmodel="${fmodel}/basidio-dmodel-plusminus-C-10-d-25-70.hd5"
amodel="${fmodel}/basidio-amodel-plusminus-C-10-d-25-70-subsampled.hd5"
imodel="${fmodel}/intron-model-C-7-d-6.hd5"
phylum="basidiomycota"

while read fungi; do
  assembly="${dataf}/Assembly/${fungi}_AssemblyScaffolds.fasta"
  introns="${dataf}/new-sequences/${fungi}-introns.fasta"
  echo "${fungi} with assembly ${assembly} and introns ${introns}"

  qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model=${dmodel},splice_site_acceptor_model=${amodel},intron_model=${imodel},phylum=${phylum},intron_source="${introns}",fungi_name="${fungi}" pipeline_perun.sh
done < "basidiomycota.txt"
