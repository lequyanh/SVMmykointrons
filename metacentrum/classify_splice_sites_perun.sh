#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=16:mem=8gb:scratch_local=2gb -v MODEL=path/to/model,DATAFILE=ascomycota_acceptor_plus_minus_testset.csv,I_WIN=70,SUBJECT=donor,IMBALANCE_RAT=0.025,CPU=16 classify_splice_sites_perun.sh

ROOT="/storage/praha1/home/$LOGNAME/"
PYTHONDIR="/storage/praha1/home/$LOGNAME/miniconda3/envs/mykointron/bin"

echo "Processing..."
trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp "${ROOT}/data/${SUBJECT}/${DATAFILE}" "${SCRATCHDIR}" || exit 1
echo "Loaded data"
cp "${ROOT}/scripts/classify-splice-sites.py" "${SCRATCHDIR}" || exit 1
cp "${ROOT}/scripts/tools.py" "${SCRATCHDIR}" || exit 1
echo "Loaded scripts"

cd "${SCRATCHDIR}" || exit 2

echo "Using model at ${MODEL}"

mkdir -p "${ROOT}/results/"

$PYTHONDIR/python classify-splice-sites.py "${DATAFILE}" "${MODEL}" "${I_WIN}" 70 "${SUBJECT}" -r "${IMBALANCE_RAT}" -c "${CPU}"

mv classify-splice-sites*.log "${ROOT}/results/"

rm -rf $SCRATCHDIR/*
