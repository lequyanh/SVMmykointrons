#!/bin/bash

DIVISION="basidiomycota"

donor_model=$1
acceptor_model=$2
train_test=$3

ROOT="/storage/praha1/home/${LOGNAME}/"

donor_pos="${ROOT}/results/${DIVISION}/donor_positions/${train_test}"
acceptor_pos="${ROOT}/results/${DIVISION}/acceptor_positions/${train_test}"

mkdir -p "${donor_pos}"
mkdir -p "${acceptor_pos}"

while read shroom_name; do
  donor_file="${ROOT}/data/${DIVISION}/donor_candidates/${train_test}/${shroom_name}_donor_cands"
  acceptor_file="${ROOT}/data/${DIVISION}/acceptor_candidates/${train_test}/${shroom_name}_acceptor_cands"

  echo "Classifying splice site candidates from ${donor_file} and ${acceptor_file}"

  # prepare files for the donor and acceptor classification results
  donor_result="${donor_pos}/${shroom_name}_results"
  acceptor_result="${acceptor_pos}/${shroom_name}_results"

  if [ -f "$donor_result" ]; then
    continue
  fi
    # classify the donors and acceptors in parallel
  qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DONOR_MODEL="${donor_model}",ACCEPTOR_MODEL="${acceptor_model}",DONOR_FILE="${donor_file}",ACCEPTOR_FILE="${acceptor_file}",DONOR_RESULT="${donor_result}",ACCEPTOR_RESULT="${acceptor_result}" pipestep_extract_positive_splices.sh

done <"${ROOT}/data/${DIVISION}/intron_${train_test}_names.txt"
