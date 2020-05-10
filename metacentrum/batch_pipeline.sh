fmodel=/storage/praha1/home/lequyanh/data/basidiomycota/models
dataf=/storage/praha1/home/lequyanh/data/reduced

dmodel="${fmodel}/dmodel-basi-intragen-d-25-C-10-win-70.hd5"
amodel="${fmodel}/amodel-basi-intragen-d-25-C-10-win-70.hd5"
imodel="${fmodel}/intron-model-C-6-d-4.hd5"
phylum="basidiomycota"

while read fungi; do
  assembly="${dataf}/Assembly/${fungi}_AssemblyScaffolds.fasta"
  introns="${dataf}/new-sequences/${fungi}-introns.fasta"
  echo "${fungi} with assembly ${assembly} and introns ${introns}"

  qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=4gb:scratch_local=2gb -v assembly_filepath="${assembly}",splice_site_donor_model=${dmodel},splice_site_acceptor_model=${amodel},intron_model=${imodel},phylum=${phylum},intron_souce="${introns}",fungi_name="${fungi}" pipeline_perun.sh
done < "basidiomycota.txt"
