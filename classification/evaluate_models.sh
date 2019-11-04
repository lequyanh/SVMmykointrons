#!/bin/bash

# bash evaluate_models.sh -f ~/Desktop/results/ -r 0.025 -t 0.25 -a 0.5 -g
# bash evaluate_models.sh -f ~/Desktop/results/Argma1-train/result.txt -r 0.025 -t 0.25 -a 0.5

g_flag=""
p_real=1
p_test=1
alpha=0

function print_usage() {
    echo "Use -g (optional) if evaluating folder of results. Outputs top metrics and the corresponding model name"
    echo "Use -f <folder/file> is either a folder of multiple results or a single result.txt file"
    echo "Use -r <real prior of minor class> e.g. 0.025"
    echo "Use -t <test prior of minor class> e.g. 0.25"
    echo "Use -a <recall/precision weight> to prioritize either recall or precision. Hihger value will favor precision"
}

# Folder with result files. Each folder having a result.txt file
while getopts 'ga:f:r:t:' flag; do
  case $flag in
    a) alpha=$OPTARG;;
    g) g_flag="true";;
    f) FOLDERS=$OPTARG ;;
    r) p_real=$OPTARG;;
    t) p_test=$OPTARG;;
    *) print_usage
       exit 1 ;;
  esac
done

function get_metrics() {
  ACC=$(cat $1 | grep  'Accuracy' | cut -d ' ' -f2)
  TPs=$(cat $1 | grep  'TP' | cut -d ' ' -f3)
  FPs=$(cat $1 | grep  'FP' | cut -d ' ' -f3)
  FNs=$(cat $1 | grep  'FN' | cut -d ' ' -f3)

  PREC=$(awk "BEGIN {print $TPs/($TPs + $FPs)}")
  RECALL=$(awk "BEGIN {print $TPs/($TPs + $FNs)}")

  echo "Metrics for model $1: TPs: $TPs, FPs: $FPs, FNs = $FNs"

  ADJ_PREC=$(awk "BEGIN {print ($p_real/$p_test)*$TPs/(($p_real/$p_test)*$TPs + (1-$p_real/1-$p_test)*$FPs)}")
  SCORE=$(awk "BEGIN {print $alpha*$ADJ_PREC + (1-$alpha)*$RECALL}")

  echo "$SCORE"
}

function best_model_metrics(){
  MAXRC=0
  MAXPREC=0
  MAXSCORE=0

  for f in $1/*.txt
  do

    get_metrics "$f"

    if (( $(echo "$RECALL > $MAXRC" |bc -l) )); then
        MAXRC=$RECALL
        MAXRC_SOURCE=$f
    fi

    if (( $(echo "$ADJ_PREC > $MAXPREC" |bc -l) )); then
        MAXPREC=$ADJ_PREC
        MAXPREC_SOURCE=$f
    fi

    if (( $(echo "$SCORE > $MAXSCORE" |bc -l) )); then
        MAXSCORE=$SCORE
        MAXSCORE_SOURCE=$f
    fi
  done
}

if [ "$g_flag" = "true" ]; then
  best_model_metrics "$FOLDERS"
  echo "Max score of $MAXSCORE at $MAXSCORE_SOURCE"
  echo "Max adjusted precision of $MAXPREC at $MAXPREC_SOURCE"
  echo "Max recall of $MAXRC at $MAXRC_SOURCE"
else
  get_metrics "$FOLDERS"
  echo "Score $SCORE, recall $RECALL, adjusted precision $ADJ_PREC"
fi



