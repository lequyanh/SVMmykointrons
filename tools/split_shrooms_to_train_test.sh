#!/bin/bash

SHROOM_NAMES_ALL=$1
NO_SHROOMS=$(wc -l < ${SHROOM_NAMES_ALL}) 
TRAIN_SIZE=$(( ${NO_SHROOMS} / 2 * 6 / 10 ))
TEST_SIZE=$(( ${NO_SHROOMS} / 2 * 4 / 10 ))

echo Number of mushrooms ${NO_SHROOMS}. Splitting to 2 train lists with ${TRAIN_SIZE} and 2 test lists with ${TEST_SIZE} classes

{
    head -n ${TRAIN_SIZE} >  splice_sites_train_names.txt
    head -n ${TEST_SIZE} > splice_sites_test_names.txt
    head -n ${TRAIN_SIZE} > intron_train_names.txt
    cat >  intron_test_names.txt
} < ${SHROOM_NAMES_ALL}
