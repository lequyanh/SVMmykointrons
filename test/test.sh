#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd ../pipeline || exit

##########
# TEST 1 #
##########
echo '===== Test simple sharding with duplicated sequences ===='
TEST_PROJECT="$SCRIPT_DIR/projects/project_test_shard_min"

bash shard_assembly.sh -p "$TEST_PROJECT" -n 2

if [ -f "$TEST_PROJECT/contigA_no_duplicates.fa" ]; then
  no_seqs=$(grep -c ">" <"$TEST_PROJECT/contigA_no_duplicates.fa")
  echo "TEST: Check correct number of sequences:"
  [[ $no_seqs -eq 4 ]] && echo True || echo False
fi

# Test the assembly got sharded into correct number of sequences
if [ -f "$TEST_PROJECT/assembly_shards/contigA_0.fa" ] \
 && [ -f "$TEST_PROJECT/assembly_shards/contigA_2.fa" ] ;then
   echo "TEST: Both shard files exist."
else
   echo "TEST: Shard files dont exist"
   exit 1
fi

# Test cleanup
rm -rf "${TEST_PROJECT:?}/assembly_shards/"
rm "${TEST_PROJECT}"/*_no_duplicates.fa

##########
# TEST 2 #
##########
echo '===== Test simple sharding and processing with duplicated sequences (neural nets models)===='
TEST_PROJECT="$SCRIPT_DIR/projects/project_test_shard_min"

bash shard_assembly.sh -p "$TEST_PROJECT" -n 2
bash batch_pipeline.sh -m nn100 -p "$TEST_PROJECT" -s plus
bash combine_results.sh -p "$TEST_PROJECT"

# Test the assembly got sharded into correct number of sequences
if [ -f "$TEST_PROJECT/results/contigA_0.fa_results_plus/pruned-contigA_0.fa" ] \
 && [ -f "$TEST_PROJECT/results/contigA_2.fa_results_plus/pruned-contigA_2.fa" ] ;then
   echo "TEST: Both results files exist."
else
   echo "TEST: Results files not present"
   exit 1
fi

# Test the assembly got sharded into correct number of sequences
if [ -f "$TEST_PROJECT/pruned_contigA_no_duplicates_plus.fa" ] ;then
   echo "TEST: Final pruned assembly file exists."
else
   echo "TEST: Final pruned assembly file not present"
   exit 1
fi

# Test cleanup
rm -rf "${TEST_PROJECT:?}/assembly_shards/"
rm -rf "${TEST_PROJECT:?}/scratchdir/"
rm -rf "${TEST_PROJECT:?}/results/"
rm "${TEST_PROJECT}"/*_no_duplicates.fa
rm "${TEST_PROJECT}"/assembly_shards.txt
rm "${TEST_PROJECT}"/pruned_*.fa
rm "${TEST_PROJECT}"/cut-coords*.csv

###########
## TEST 3 #
###########
echo '===== Test simple sharding and processing with duplicated sequences (neural nets models)===='
TEST_PROJECT="$SCRIPT_DIR/projects/project_test_shard_min"

bash shard_assembly.sh -p "$TEST_PROJECT" -n 1
bash batch_pipeline.sh -m nn100 -p "$TEST_PROJECT" -s both
bash combine_results.sh -p "$TEST_PROJECT"

# Test the assembly got sharded into correct number of sequences
if [ -f "$TEST_PROJECT/results/contigA_0.fa_results_plus/pruned-contigA_0.fa" ] \
 && [ -f "$TEST_PROJECT/results/contigA_3.fa_results_plus/pruned-contigA_3.fa" ] ;then
   echo "TEST: Both results files exist."
else
   echo "TEST: Results files not present"
   exit 1
fi

# Test the assembly got sharded into correct number of sequences
if [ -f "$TEST_PROJECT/pruned_contigA_no_duplicates_plus.fa" ] \
 && [ -f "$TEST_PROJECT/pruned_contigA_no_duplicates_minus.fa" ];then
   echo "TEST: Final pruned assembly file exists."
else
   echo "TEST: Final pruned assembly file not present"
   exit 1
fi

# Test cleanup
rm -rf "${TEST_PROJECT:?}/assembly_shards/"
rm -rf "${TEST_PROJECT:?}/scratchdir/"
rm -rf "${TEST_PROJECT:?}/results/"
rm "${TEST_PROJECT}"/*_no_duplicates.fa
rm "${TEST_PROJECT}"/assembly_shards.txt
rm "${TEST_PROJECT}"/pruned_*.fa
rm "${TEST_PROJECT}"/cut-coords*.csv

echo "SUCCESS"
