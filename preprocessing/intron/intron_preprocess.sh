#!/bin/bash

# System settings:
#  - path to python
PYTHON=~/anaconda3/envs/mykointron/bin/python
#  - number of CPUs to be used in total
#    the minimal value is 2, although 16 or more is recommended
NUMBER_CPUS=10
#  - memory limit for caching (in MB)
#    this is used in the intron prediction
CACHE_LIMIT=1024
# -------------------------------------------------------------

# Pipeline settings:
#  - splice site dimers
DONOR="GT"
ACCEPTOR="AG"
#  - window size used when extracting splice sites sequences
#    the size is equal to the size of window used to train the models
DONOR_LWINDOW=70
DONOR_RWINDOW=70
ACCEPTOR_LWINDOW=70
ACCEPTOR_RWINDOW=70

# Number of scaffolds taken from each shroom for intron model training
NO_SCAFF=5
#  - range of intron lengths
#    considered when extracting introns from the positions of the positively classified splice sites
INTRON_MIN_LENGTH=10
INTRON_MAX_LENGTH=600
#  - order of the spectrum kernel
#    it is used in the intron prediction and it must be equal to the order used while training
SPECT_KERNEL_ORDER=2
#  - name of the file that contains positions of alleged introns
#    it will be created after the splice sites classification
INTRON_POSITIONS_FILE="intron-positions-dataset.csv"
#  - name of the file that contains extracted intron sequences
INTRON_FILE="intron-dataset.csv"
#  - name of the file for intron classification results
INTRON_RESULT="intron-result.csv"

# -------------------------------------------------------------

# Derived variables:
# regex used to determine splice site sequences
donor_regex=";[ACGT]{$DONOR_LWINDOW}$DONOR[ACGT]{$DONOR_RWINDOW}$"
acceptor_regex=";[ACGT]{$ACCEPTOR_LWINDOW}$ACCEPTOR[ACGT]{$ACCEPTOR_RWINDOW}$"
# regex to filter positively classified splice sites
positive_splice_sites=";1$"
# regex to filter positively classified introns
positive_introns=";1$"
# -------------------------------------------------------------

# Pipeline inputs:
#  - Assembly file (FASTA)
assebmlies_loc=$1   #/home/anhvu/Desktop/mykointrons-data/data/Assembly
#  - donor splice site prediction model
splice_site_donor_model=$2
#  - acceptor splice site prediction model
splice_site_acceptor_model=$3
#  - intron prediction model
intron_model=$4

if [ $# -ne 4 ]
then
    echo "Arguments expected: ASSEMBLY DONOR_MODEL ACCEPTOR_MODEL INTRON_MODEL"
    exit 1
fi

# -------------------------------------------------------------

# The process (roughly):
#  1. Find all donors (GT) and use them to create a donor dataset.
#  2. Find all acceptors (AG) and use them to create an acceptor dataset.
#  3. Run two classification tasks in parallel:
#     a) classify donors (given the donor dataset),
#     b) classify acceptors (given the acceptor dataset).
#  4. Wait until both the classification tasks finish. Then find all possible pairs GT-AG given
#     positively classified donors and acceptors, and use them to create an intron dataset.
#  5. Run a classification task to identify introns.
# -------------------------------------------------------------

# exit after any error
set -e
# enable pipe fail
set -o pipefail

#function generate_splice_site_candidates() {
#  train_test=$1
#
#  donors_cands_loc="donor_candidates/${train_test}"
#  accepts_cands_loc="acceptor_candidates/${train_test}"
#
#  mkdir -p "${donors_cands_loc}"
#  mkdir -p "${accepts_cands_loc}"
#
#  while read shroom_name; do
#  echo "$shroom_name"
#
#  assembly_filepath="${assebmlies_loc}/${shroom_name}_AssemblyScaffolds.fasta"
#  echo "Extracting donors and acceptors from [$assembly_filepath]..."
#
#  donor_file="${donors_cands_loc}/${shroom_name}_donor_cands"
#  acceptor_file="${accepts_cands_loc}/${shroom_name}_acceptor_cands"
#
#  # prepare files for the donor and acceptor datasets
#  echo "scaffold;position;sequence" > "$donor_file"
#  echo "scaffold;position;sequence" > "$acceptor_file"
#
#  # extract the splice site sequences and then use awk to separate the sequences into
#  # the donor and acceptor dataset, respectively
#  $PYTHON extract-donor-acceptor.py $assembly_filepath \
#                                    $DONOR $ACCEPTOR \
#                                    $DONOR_LWINDOW $DONOR_RWINDOW \
#                                    $ACCEPTOR_LWINDOW $ACCEPTOR_RWINDOW \
#                                    $NO_SCAFF \
#          | gawk -v donor="$donor_regex" \
#                -v acceptor="$acceptor_regex" \
#                -v donor_file=$donor_file \
#                -v acceptor_file=$acceptor_file \
#                '$0 ~ donor {print >> donor_file} $0 ~ acceptor {print >> acceptor_file}'
#
#  echo "Donors extracted to [$donor_file]."
#  echo "Acceptors extracted to [$acceptor_file]."
#  echo ""
#
#done < ../shroom_split/"intron_${train_test}_names.txt"
#}
#
#echo "Generate acceptor and donor candidates for intron training"
#generate_splice_site_candidates "train"
#echo "Generate acceptor and donor candidates for intron testing"
#generate_splice_site_candidates "test"

# determine the number of CPUs available for each classification task
donor_cpus=$((NUMBER_CPUS/2))
acceptor_cpus=$((NUMBER_CPUS-donor_cpus))

echo "Starting classification of splice sites with [$donor_cpus/$acceptor_cpus] CPUs..."

function get_positive_splice_site_candidates() {
  train_test=$1

  donor_pos="donor_positions/${train_test}"
  acceptor_pos="acceptor_positions/${train_test}"

  mkdir -p "${donor_pos}"
  mkdir -p "${acceptor_pos}"

  while read shroom_name; do
    donor_file="donor_candidates/${train_test}/${shroom_name}_donor_cands"
    acceptor_file="acceptor_candidates/${train_test}/${shroom_name}_acceptor_cands"

    echo "Classifying splice site candidates from ${donor_file} and ${acceptor_file}"

    # prepare files for the donor and acceptor classification results
    donor_result="${donor_pos}/${shroom_name}_results_txt"
    acceptor_result="${acceptor_pos}/${shroom_name}_results.txt"

    echo "scaffold;position" > "$donor_result"
    echo "scaffold;position" > "$acceptor_result"

    echo "Positive splice site positions will be saved to ${donor_result} and ${acceptor_result}"

    # classify the donors and acceptors in parallel
    # keep only the positively classified samples
    # keep only the columns `scaffold`, and `position` (1st and 2nd)
    $PYTHON classify-splice-sites.py $donor_file $splice_site_donor_model \
                                    $DONOR_RWINDOW $DONOR_LWINDOW \
                                    "donor" \
                                    -c $donor_cpus \
            | grep $positive_splice_sites \
            | cut -d ';' -f -2 >> $donor_result &
    classify_donor_pid=$!
    $PYTHON classify-splice-sites.py $acceptor_file $splice_site_acceptor_model \
                                    $ACCEPTOR_LWINDOW $ACCEPTOR_RWINDOW \
                                    "acceptor" \
                                    -c $acceptor_cpus \
            | grep $positive_splice_sites \
            | cut -d ';' -f -2 >> $acceptor_result &
    classify_acceptor_pid=$!

    # wait for both the classification tasks to finish
    wait $classify_donor_pid $classify_acceptor_pid
  done < ../shroom_split/"intron_${train_test}_names.txt"
}

echo "Classifying splice site candidates (keeping positions of positive donors/acceptors). Creating traning dataset"
get_positive_splice_site_candidates "train"
echo "Classifying splice site candidates (keeping positions of positive donors/acceptors). Creating testing dataset"
get_positive_splice_site_candidates "test"


## given the splice sites classification, output positions of possible introns
#$PYTHON generate-pairs.py $DONOR_RESULT $ACCEPTOR_RESULT \
#                          $INTRON_MIN_LENGTH $INTRON_MAX_LENGTH > $INTRON_POSITIONS_FILE
#
#                          echo "Intron positions are in [$INTRON_POSITIONS_FILE]."
#echo "Extracting introns from the positions..."
#
## extract introns from the positions from the previous step
#$PYTHON extract-introns.py $assembly_filepath $INTRON_POSITIONS_FILE > $INTRON_FILE
#
#echo "Intron sequences are extracted in [$INTRON_FILE]."
#echo ""
#
## prepare a file for the intron classification results
#echo "scaffold;start;end" > $INTRON_RESULT
#
## classify the introns
## keep only the positively classified samples
## keep only the columns `scaffold`, `start`, and `end` (1st, 2nd, and 3rd)
#$PYTHON classify-introns.py -c $NUMBER_CPUS $INTRON_FILE $intron_model $SPECT_KERNEL_ORDER \
#        | grep $positive_introns \
#        | cut -d ';' -f -3 >> $INTRON_RESULT
#
#echo "Detected introns are in [$INTRON_RESULT]."
#
