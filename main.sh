# Path to project folder containing an fasta assembly to purge
project_path=$1
# Type of model to cut the introns with (determines other pipeline settings)
model_settings=$2 #nn100 of svmb

assembly_path=$(find "$project_path" -maxdepth 1 -name "*.fasta" -o -name "*.fa")

# Sharding settings (metagenom can be split to chunks for parallel processing)
NUM_FILES=300


##############################
# SHARDING THE MAIN ASSEMBLY #
##############################
cd ./pipeline || exit

total_lines=$(wc -l <"${assembly_path}")
((lines_per_file = (total_lines + NUM_FILES - 1) / NUM_FILES))

if ((lines_per_file % 2 != 0)); then
  ((lines_per_file = lines_per_file + 1))
fi

# Split the actual file, maintaining lines.
assembly_name=$(basename assembly_path) #TODO remove extension
split --lines="${lines_per_file}" "${assembly_path}" "${assembly_name}" -d -a 3

# Debug information
echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"

# Move metagenom shards to a separate folder
ASSEMBL_SHARDS_DIR='assembly_shards'
mkdir -p $ASSEMBL_SHARDS_DIR
find . -maxdepth 1 -type f -name "${assembly_name}[0-9]*.fa" -type f -exec mv -f {} $ASSEMBL_SHARDS_DIR \;

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

bash batch_pipeline.sh -m "$model_settings" -l "$fasta_list_file" -d $ASSEMBL_SHARDS_DIR -s plus
mkdir -p ${forw_strand_result_loc}
mv assembly_shard*/ ${forw_strand_result_loc}

bash batch_pipeline.sh -m "$model_settings" -l "$fasta_list_file" -d $ASSEMBL_SHARDS_DIR -s minus
mkdir -p ${rev_strand_result_loc}
mv assembly_shard*/ ${rev_strand_result_loc}

merge_cuts ${forw_strand_result_loc} $OUT_FORW
merge_cuts ${rev_strand_result_loc} $OUT_REV
