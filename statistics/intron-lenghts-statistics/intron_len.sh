#!/bin/bash

#aggreg_filename="$(pwd)/ascomycota-intron-lens.txt"
#new_sequences_loc="/home/anhvu/Desktop/mykointrons-data/new-sequences/"
#
#for f in ${new_sequences_loc}/*/; do
#  cd "$f" || exit
#  cat "$(basename ${f})-introns-len.txt" > $aggreg_filename
#done
PYTHON=~/anaconda3/envs/mykointron/bin/python

INTRON_MIN_LENGTH=10
INTRON_MAX_LENGTH=200
BASEDIR="/home/anhvu/PycharmProjects/mycointrons"
FASTA_DIR="/home/anhvu/Desktop/mykointrons-data/reduced/"


function extract_and_label() {
  assembly_filepath=$1
  intron_source=$2
  donor_result=$3
  acceptor_result=$4
  strand=$5
  INTRON_FILE=$6

  INTRON_POSITIONS_FILE='temp_intron_pos_file.csv'

  # given the splice sites classification, output positions of possible introns
  $PYTHON generate-pairs.py $donor_result $acceptor_result \
    $INTRON_MIN_LENGTH $INTRON_MAX_LENGTH >$INTRON_POSITIONS_FILE

  echo "Intron positions are in [$INTRON_POSITIONS_FILE]."
  echo "Extracting introns from the positions..."

  # extract introns from the positions from the previous step
  $PYTHON extract-introns.py "$assembly_filepath" $INTRON_POSITIONS_FILE "$strand" >$INTRON_FILE

  echo "Labeling intron dataset"
  $PYTHON label_introns.py "$INTRON_FILE" "$intron_source" $INTRON_MIN_LENGTH $INTRON_MAX_LENGTH "$strand"

  rm $INTRON_POSITIONS_FILE
}

fungi_list=$1
strand=$2

while read p
do
    fungi=$(echo $p | sed "s/_results//")
    echo $fungi

    assembly="${FASTA_DIR}/Assembly/${fungi}_AssemblyScaffolds.fasta"
    introns="${FASTA_DIR}/new-sequences/${fungi}-introns.fasta"

    donor_results="${BASEDIR}/intragen/sample_results_basidio_revisited/${fungi}_results/splice-site-donor-result.csv"
    acceptor_results="${BASEDIR}/intragen/sample_results_basidio_revisited/${fungi}_results/splice-site-acceptor-result.csv"

    out_introns_dataset="${BASEDIR}/statistics/intron-lenghts-statistics/basi_introns_dataset/${fungi}_introns.csv"

    extract_and_label $assembly $introns $donor_results $acceptor_results $strand $out_introns_dataset
done < ${fungi_list}