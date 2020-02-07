#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DONOR_MODEL=dmodel.hd5,ACCEPTOR_MODEL=amodel.hd5,DONOR_FILE=train_set.csv,ACCEPTOR_FILE=train_set.csv,DONOR_RESULT=dresult.csv,ACCEPTOR_RESULT=aresult.csv pipe_extract_positive_splices.sh
DONOR_LWINDOW=70
DONOR_RWINDOW=70
ACCEPTOR_LWINDOW=70
ACCEPTOR_RWINDOW=70

positive_splice_sites=";1$"

ROOT="/storage/praha1/home/${LOGNAME}/"
PYTHONDIR="/storage/praha1/home/${LOGNAME}/miniconda3/envs/mykointron/bin"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens
cp "${ROOT}/scripts/classify-splice-sites.py" "${SCRATCHDIR}" || exit 1
cp "${ROOT}/scripts/tools.py" "${SCRATCHDIR}" || exit 1
echo "Loaded scripts"

cd "${SCRATCHDIR}" || exit 2

echo "scaffold;position" >"${DONOR_RESULT}"
echo "scaffold;position" >"${ACCEPTOR_RESULT}"

echo "Positive splice site positions will be saved to ${DONOR_RESULT} and ${ACCEPTOR_RESULT}"

# classify the donors and acceptors in parallel
# keep only the positively classified samples
# keep only the columns `scaffold`, and `position` (1st and 2nd)
$PYTHONDIR/python classify-splice-sites.py "${DONOR_FILE}" "${DONOR_MODEL}" \
  $DONOR_RWINDOW $DONOR_LWINDOW \
  "donor" \
  -c 5 |
  grep $positive_splice_sites |
  cut -d ';' -f -2 >>"${DONOR_RESULT}" &
classify_donor_pid=$!
$PYTHONDIR/python classify-splice-sites.py "${ACCEPTOR_FILE}" "${ACCEPTOR_MODEL}" \
  "${ACCEPTOR_LWINDOW}" "${ACCEPTOR_RWINDOW}" \
  "acceptor" \
  -c 5 |
  grep $positive_splice_sites |
  cut -d ';' -f -2 >>"${DONOR_RESULT}" &
classify_acceptor_pid=$!

# wait for both the classification tasks to finish
wait $classify_donor_pid $classify_acceptor_pid
echo "Done"

rm -rf $SCRATCHDIR/*
