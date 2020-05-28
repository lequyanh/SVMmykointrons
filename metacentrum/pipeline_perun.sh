#!/bin/bash

# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v assembly_filepath=#Assembly,splice_site_donor_model=#Dmodel,splice_site_acceptor_model=#Amodel,intron_model=#Imodel,phylum=#phylum,intron_source=#introns-fasta,fungi_name=#fungi pipeline_perun.sh

# System settings:
#  - path to python
PYTHON="/storage/praha1/home/${LOGNAME}/miniconda3/envs/mykointron/bin/python"
#  - number of CPUs to be used in total
#    the minimal value is 2, although 16 or more is recommended
NUMBER_CPUS=16

ROOT="/storage/praha1/home/${LOGNAME}"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens
cp "${ROOT}"/scripts/pipeline/* "${SCRATCHDIR}" || exit 1
echo "Loaded scripts"

cd "${SCRATCHDIR}" || exit 2

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
#  - range of intron lengths
#    considered when extracting introns from the positions of the positively classified splice sites
STRAND="-"
INTRON_MIN_LENGTH=40
INTRON_MAX_LENGTH=100
#  - order of the spectrum kernel
#    it is used in the intron prediction and it must be equal to the order used while training
SPECT_KERNEL_ORDER_BASI=4
SPECT_KERNEL_ORDER_ASCO=6
#  - files with intron lengths to build a probability distribution over intron lengths
#    used as the cutting step to decide which candidate to cut in case of overlap
INTRON_LENGTH_DATA_BASI="basidiomycota-intron-lens.txt"
INTRON_LENGTH_DATA_ASCO="ascomycota-intron-lens.txt"

# OUTPUT FILES
#  - names of the splice site datasets
#    it will be created before the splice sites classification
DONOR_FILE="splice-site-donor-dataset.csv"
ACCEPTOR_FILE="splice-site-acceptor-dataset.csv"
#  - names of the files for splice site classification results
DONOR_RESULT="splice-site-donor-result.csv"
ACCEPTOR_RESULT="splice-site-acceptor-result.csv"
#  - name of the file that contains positions of alleged introns
#    it will be created after the splice sites classification
INTRON_POSITIONS_FILE="intron-positions-dataset.csv"
#  - name of the file that contains extracted intron sequences
INTRON_FILE="intron-dataset.csv"
#  - name of the file for intron classification results
INTRON_RESULT="intron-result.csv"
CUT_COORDS_FILE="cut-coords.csv"
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

## Pipeline inputs:
##  - Assembly file (FASTA)
#assembly_filepath=$1 # ~/Desktop/mykointrons-data/data/Assembly/Kocim1_AssemblyScaffolds.fasta
##  - donor splice site prediction model
#splice_site_donor_model=$2
##  - acceptor splice site prediction model
#splice_site_acceptor_model=$3
##  - intron prediction model
#intron_model=$4
##  - basiomycota or ascomycota. For adjusting parameters
#phylum=$5
##  - in case of validation
#intron_source=$6 #"/home/anhvu/Desktop/mykointrons-data/new-sequences/Thega1/Thega1-introns.fasta"
#
#fungi_name=$7

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

function init() {
  case $phylum in
  "ascomycota")
    echo "Loading settings for ${phylum}"
    spect_kernel_order=${SPECT_KERNEL_ORDER_ASCO}
    intron_lens_data=${INTRON_LENGTH_DATA_ASCO}
    ;;

  "basidiomycota")
    echo "Loading settings for ${phylum}"
    spect_kernel_order=${SPECT_KERNEL_ORDER_BASI}
    intron_lens_data=${INTRON_LENGTH_DATA_BASI}
    ;;

  *)
    echo "Phylum not known, exiting" &
    exit
    ;;
  esac
}

function extract_donor_acceptor_step() {
  echo "Extracting donors and acceptors from [$assembly_filepath] on $STRAND strand..."

  # prepare files for the donor and acceptor datasets
  echo "scaffold;position;sequence" > $DONOR_FILE
  echo "scaffold;position;sequence" > $ACCEPTOR_FILE

  $PYTHON extract-donor-acceptor.py "${assembly_filepath}" \
    $DONOR $ACCEPTOR \
    $DONOR_LWINDOW $DONOR_RWINDOW \
    $ACCEPTOR_LWINDOW $ACCEPTOR_RWINDOW \
    $STRAND \
|gawk -v donor="${donor_regex}" \
      -v acceptor="${acceptor_regex}" \
      -v donor_file=$DONOR_FILE \
      -v acceptor_file=$ACCEPTOR_FILE \
      '$0 ~ donor {print >> donor_file} $0 ~ acceptor {print >> acceptor_file}'

  echo "Donors extracted to [$DONOR_FILE]."
  echo "Acceptors extracted to [$ACCEPTOR_FILE]."
  echo ""
}

function classify_splice_sites_step() {
  # prepare files for the donor and acceptor classification results
  echo "scaffold;position" > $DONOR_RESULT
  echo "scaffold;position" > $ACCEPTOR_RESULT

  echo "Positive splice site positions will be saved to ${DONOR_RESULT} and ${ACCEPTOR_RESULT}"

  donor_cpus=$((NUMBER_CPUS / 2))
  acceptor_cpus=$((NUMBER_CPUS - donor_cpus))

  # classify the donors and acceptors in parallel
  # keep only the positively classified samples
  # keep only the columns `scaffold`, and `position` (1st and 2nd)
  echo "Starting classification of splice sites with [$donor_cpus/$acceptor_cpus] CPUs..."

  $PYTHON classify-splice-sites.py $DONOR_FILE "$splice_site_donor_model" \
    $DONOR_RWINDOW $DONOR_LWINDOW \
    "donor" -c $donor_cpus |
    grep $positive_splice_sites |
    cut -d ';' -f -2 >>$DONOR_RESULT &
  classify_donor_pid=$!

  $PYTHON classify-splice-sites.py $ACCEPTOR_FILE "$splice_site_acceptor_model" \
    $ACCEPTOR_LWINDOW $ACCEPTOR_RWINDOW \
    "acceptor" -c $acceptor_cpus |
    grep $positive_splice_sites |
    cut -d ';' -f -2 >>$ACCEPTOR_RESULT &
  classify_acceptor_pid=$!

  # wait for both the classification tasks to finish
  wait $classify_donor_pid $classify_acceptor_pid
  echo "Done"
}

function pair_splice_sites_and_extract_introns_step() {
  echo "Getting intron positions..."

  # given the splice sites classification, output positions of possible introns
  $PYTHON generate-pairs.py $DONOR_RESULT $ACCEPTOR_RESULT \
    $INTRON_MIN_LENGTH $INTRON_MAX_LENGTH >$INTRON_POSITIONS_FILE

  echo "Intron positions are in [$INTRON_POSITIONS_FILE]."
  echo "Extracting introns from the positions..."

  # extract introns from the positions from the previous step
  $PYTHON extract-introns.py "$assembly_filepath" $INTRON_POSITIONS_FILE $STRAND > $INTRON_FILE

  echo "Intron sequences are extracted in [$INTRON_FILE]."
  echo ""

  echo "First 20 most frequent intron candidates:"
#  cat $INTRON_FILE | cut -d ';' -f4 | sort | uniq -c | grep -v '1 ' | sort -r | head -20

}

function classify_introns_step() {
  echo "Starting classification of the introns using spectral kernel order ${spect_kernel_order}"

  # prepare a file for the intron classification results
  echo "scaffold;start;end" > $INTRON_RESULT

  # classify the introns
  # keep only the positively classified samples
  # keep only the columns `scaffold`, `start`, and `end` (1st, 2nd, and 3rd)
  $PYTHON classify-introns.py -c $NUMBER_CPUS $INTRON_FILE "$intron_model" $spect_kernel_order |
    grep $positive_introns |
    cut -d ';' -f -3 >> $INTRON_RESULT

  echo "Detected introns are in [$INTRON_RESULT]."
}

function validate_introns_step() {
  echo "Labeling intron dataset"
  $PYTHON label_introns.py "$INTRON_FILE" "$intron_source" $INTRON_MIN_LENGTH $INTRON_MAX_LENGTH $STRAND

  echo "Starting validation of the intron dataset using spectral kernel order ${spect_kernel_order}"
  $PYTHON classify-introns.py -c $NUMBER_CPUS $INTRON_FILE "$intron_model" $spect_kernel_order |
    cut --complement -d ';' -f4 >> $INTRON_RESULT
}

function cut_introns_step() {
  echo "Cutting introns. Intron length distribution derived from file ${intron_lens_data}"

  echo "scaffold;start;end" > ${CUT_COORDS_FILE}
  $PYTHON prune_probabilistic.py "${assembly_filepath}" ${INTRON_RESULT} ${intron_lens_data} $STRAND >> "${CUT_COORDS_FILE}"
}

init
extract_donor_acceptor_step
classify_splice_sites_step
pair_splice_sites_and_extract_introns_step

## ===================================== VALIDATE =====================================
validate_introns_step
## =================================== or CLASIFY =====================================
#classify_introns_step

cut_introns_step

echo "Pipeline log for file ${assembly_filepath}" > pipeline.output
for p in ./*.log; do
  {
    echo "================================== ${p} =========================================="
    cat "$p"
  } >> pipeline.output
done

result_dir="${fungi_name}_results"
mkdir "${result_dir}"

mv "${CUT_COORDS_FILE}" "${result_dir}/"
mv "${ACCEPTOR_RESULT}" "${result_dir}/"
mv "${DONOR_RESULT}" "${result_dir}/"
mv "${INTRON_RESULT}" "${result_dir}/"

zip -r "${result_dir}.zip" "${result_dir}/"
mv "${result_dir}.zip" "${ROOT}/results/"

rm -rf $SCRATCHDIR/*