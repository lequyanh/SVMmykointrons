#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v ORDER=4,C_POS=3,INPUT=train_set.csv,CPU=10,TEST=0.2 train-introns.sh

ROOT="/storage/praha1/home/${LOGNAME}/"
PYTHONDIR="/storage/praha1/home/${LOGNAME}/miniconda3/envs/mykointron/bin"

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp "${ROOT}/data/introns/${INPUT}" "${SCRATCHDIR}" || exit 1
echo "Data loaded"

cp "${ROOT}/scripts/train-introns.py" "${SCRATCHDIR}" || exit 1
echo "Libraries loaded"

cd "${SCRATCHDIR}" || exit 2

cat *.csv > input.csv

$PYTHONDIR/python train-introns.py input.csv "${ORDER}" "${C_POS}" -t "${TEST}" -c "${CPU}"

mkdir -p "${ROOT}/results"

RESULTDIR="train-introns-C-${C_POS}-d-${ORDER}-${PBS_JOBID}"
mkdir "${RESULTDIR}"

mv train-introns.*.log "${RESULTDIR}/"
mv model.hd5 "${RESULTDIR}/"

tar -czvf "${RESULTDIR}.tar.gz" "${RESULTDIR}"

cp "${RESULTDIR}.tar.gz" "${ROOT}/results/"

rm -rf $SCRATCHDIR/*
