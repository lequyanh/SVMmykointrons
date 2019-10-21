#!/bin/bash

SHROOM_NAMES_ALL=$1

INT_TRAIN_NAMES=$2
INT_TEST_NAMES=$3
SPLICE_TRAIN_NAMES=$4
SPLICE_TEST_NAMES=$5

NO_SHROOMS=$(wc -l < ${SHROOM_NAMES_ALL}) 
TRAIN_SIZE=$(( ${NO_SHROOMS} / 2 * 6 / 10 ))
TEST_SIZE=$(( ${NO_SHROOMS} / 2 * 4 / 10 ))

echo Number of mushrooms ${NO_SHROOMS}. Splitting to 2 train lists with ${TRAIN_SIZE} and 2 test lists with ${TEST_SIZE} classes

{
    head -n ${TRAIN_SIZE} >  $SPLICE_TRAIN_NAMES
    head -n ${TEST_SIZE} > $SPLICE_TEST_NAMES
    head -n ${TRAIN_SIZE} > $INT_TRAIN_NAMES
    cat >  $INT_TEST_NAMES
} < ${SHROOM_NAMES_ALL}
