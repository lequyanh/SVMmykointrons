#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v MODEL_FOLDER=path/to/models,DATAFILE=Mycrab_acceptors.csv,I_WIN=70,SUBJECT=donor,IMBALANCE_RAT=0.025,CPU=10 classify-splice-sites-cloud.sh

ROOT="/storage/praha1/home/$LOGNAME/"
PYTHONDIR="/storage/praha1/home/$LOGNAME/miniconda3/envs/mykointron/bin"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp "${ROOT}/data/${SUBJECT}/${DATAFILE}" "${SCRATCHDIR}" || exit 1
echo "Loaded data"
cp "${ROOT}/scripts/classify-splice-sites.py" "${SCRATCHDIR}" || exit 1
cp "${ROOT}/scripts/tools.py" "${SCRATCHDIR}" || exit 1
echo "Loaded scripts"

cd "${SCRATCHDIR}" || exit 2

MODEL="${MODEL_FOLDER}/model.hd5"

echo "Using model at ${MODEL}"

$PYTHONDIR/python classify-splice-sites.py "${DATAFILE}" "${MODEL}" "${I_WIN}" 70 "${SUBJECT}" -r "${IMBALANCE_RAT}" -c "${CPU}"

mkdir -p "${ROOT}/results/"

mv classify-splice-sites*.log "${MODEL_FOLDER}"

rm -rf $SCRATCHDIR/*
