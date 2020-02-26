#!/bin/bash

# how to run
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v PHYLUM=ascomycota,DATAFILE=data.csv,LWINDOW=60,RWINDOW=80,SITE=donor,DEGREE=30,C=1,CPU=10 train-splice-site-meta.sh

ROOT="/storage/praha1/home/$LOGNAME/"
PYTHONDIR="/storage/praha1/home/$LOGNAME/miniconda3/envs/mykointron/bin"

ASCO_IMBALANCE=0.0033
BASIDIO_IMBALANCE=0.022

trap 'clean_scratch' TERM EXIT

# $SCRATCHDIR = path to the workstation where the computation happens

cp "${ROOT}/data/${SITE}/${DATAFILE}" "${SCRATCHDIR}" || exit 1
echo "Loaded data"

cp "${ROOT}/scripts/train-splice-sites.py" "${SCRATCHDIR}" || exit 1
cp "${ROOT}/scripts/tools.py" "${SCRATCHDIR}" || exit 1
echo "Loaded libraries"

cd "${SCRATCHDIR}" || exit 2

cat *.csv > input.csv

model_filename="model.hd5"
[[ $PHYLUM = "ascomycota" ]] && imbal_ratio=${ASCO_IMBALANCE} || imbal_ratio=${BASIDIO_IMBALANCE}
echo "Imbalance ratio of ${imbal_ratio}"

$PYTHONDIR/python train-splice-sites.py input.csv ${model_filename} "${LWINDOW}" "${RWINDOW}" "${SITE}" "${DEGREE}" "${C}" -r ${imbal_ratio} -c "${CPU}"

mkdir -p "${ROOT}"/results

RESULTDIR="train-$SITE-$DATAFILE-C-$C-d-$DEGREE-win-$LWINDOW-$RWINDOW-$PBS_JOBID"
mkdir "${RESULTDIR}"

mv train-splice-sites-.*.log "${RESULTDIR}/"
mv model.hd5 "${RESULTDIR}/"

tar -czvf "${RESULTDIR}.tar.gz" "${RESULTDIR}"

cp "${RESULTDIR}.tar.gz" "${ROOT}/results/"

rm -rf $SCRATCHDIR/*
