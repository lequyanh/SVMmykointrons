#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v MODEL_FOLDER=path/to/models,DATAFILE=intron_test_set.csv,ORDER=4,CPU=10 validate_gridmodel_intron_meta.sh

ROOT="/storage/praha1/home/$LOGNAME/"
PYTHONDIR="/storage/praha1/home/$LOGNAME/miniconda3/envs/mykointron/bin"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp "${ROOT}/data/introns/${DATAFILE}" "${SCRATCHDIR}" || exit 1
echo "Loaded data"
cp "${ROOT}/scripts/classify-introns.py" "${SCRATCHDIR}" || exit 1
cp "${ROOT}/scripts/tools.py" "${SCRATCHDIR}" || exit 1
echo "Loaded scripts"

cd "${SCRATCHDIR}" || exit 2

MODEL="${MODEL_FOLDER}/model.hd5"

$PYTHONDIR/python classify-introns.py "${DATAFILE}" "${MODEL}" "${ORDER}" -c "${CPU}"
cat classify-introns*.log >"${MODEL_FOLDER}/validation-results.txt"

rm classify-introns*.log

mv classify-splice-sites*.log "${MODEL_FOLDER}"

rm -rf $SCRATCHDIR/*
