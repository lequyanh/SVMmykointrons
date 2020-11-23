METAGENOM_FILEPATH=$1
# Type of model to cut the introns with (determines other pipeline settings)
MODEL_SETTINGS=$2 #nn100 for instance

# Sharding settings (metagenom can be split to chunks for parallel processing)
NUM_FILES=300
# Output file names
OUT_FORW='metagenom-cut-coords.csv'
OUT_REV='metagenom-cut-coords-reverse.csv'

##########################
# SHARDING THE METAGENOM #
##########################
cd ./pipeline || exit

total_lines=$(wc -l <"${METAGENOM_FILEPATH}")
((lines_per_file = (total_lines + NUM_FILES - 1) / NUM_FILES))

if ((lines_per_file % 2 != 0)); then
  ((lines_per_file = lines_per_file + 1))
fi

# Split the actual file, maintaining lines.
split --lines=${lines_per_file} "${METAGENOM_FILEPATH}" metagenom_shard -d -a 3

# Debug information
echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"

# Move metagenom shards to a separate folder
META_SHARDS_F='meta_shards'
mkdir -p $META_SHARDS_F
find . -maxdepth 1 -type f -name 'metagenom_shard[0-9]*' -type f -exec mv -f {} $META_SHARDS_F \;

# list of metagenom shards to cut
fasta_list_file="metagenom_shards.txt"
ls $META_SHARDS_F >$fasta_list_file
# Can be replaced by custom list of shards (in case we need to re-compute some etc.)
#fasta_list_file="custom_metagenom_shard_list.txt"
#todo remove this
head -3 $fasta_list_file >"reduced-${fasta_list_file}"
fasta_list_file="reduced-${fasta_list_file}"
#todo remove this

###################
# GET INTRON CUTS #
###################
# Function for merging cuts into one file
function merge_cuts() {
  results_folder=$1
  out=$2

  echo "scaffold;start;end" >"${out}"

  for f in "${results_folder}"/*/; do
    # Read all lines except for the first one containing headers
    tail -n +2 "$f/cut-coords.csv" >>"${out}"
  done
}

# Cut both strands
forw_strand_result_loc="metagenom_forw_stand_res"
rev_strand_result_loc="metagenom_rev_stand_res"

bash batch_pipeline.sh -m "$MODEL_SETTINGS" -l "$fasta_list_file" -d $META_SHARDS_F -s + -f
mkdir -p ${forw_strand_result_loc}
mv metagenom_shard*/ ${forw_strand_result_loc}

bash batch_pipeline.sh -m "$MODEL_SETTINGS" -l "$fasta_list_file" -d $META_SHARDS_F -s - -f
mkdir -p ${rev_strand_result_loc}
mv metagenom_shard*/ ${rev_strand_result_loc}

merge_cuts ${forw_strand_result_loc} $OUT_FORW
merge_cuts ${rev_strand_result_loc} $OUT_REV

#mv $OUT_FORW ../annotation
#mv $OUT_REV ../annotation

#####################
# ANNOTATE THE CUTS #
#####################
cd ../annotation || exit
bash INTRONY_PIPELINE.sh "$METAGENOM_FILEPATH" $OUT_FORW $OUT_REV
