 #!/bin/bash
a_flag=""

function print_usage() {
    echo "Use -a (optional) if evaluating folder of results. Outputs top metrics and the corresponding model name"
    echo "Use -f <folder/file> is either a folder of multiple results or a single result.txt file"
}

# Folder with result files. Each folder having a result.txt file
while getopts 'af:' flag; do
  case $flag in
    a) a_flag="true" ;;
    f) FOLDERS=$OPTARG ;;
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

  ADJ_RECALL=0
}

function best_model_metrics(){
  MAXRC=0
  MAXPREC=0
  MAXACC=0

  for f in $1/*/
  do
    cd $f
    get_metrics result.txt

    if (( $(echo "$RECALL > $MAXRC" |bc -l) )); then
        MAXRC=$RECALL
        MAXRC_SOURCE=$f
    fi

    if (( $(echo "$PREC > $MAXPREC" |bc -l) )); then
        MAXPREC=$PREC
        MAXPREC_SOURCE=$f
    fi

    if (( $(echo "$ACC > $MAXACC" |bc -l) )); then
        MAXACC=$ACC
        MAXACC_SOURCE=$f
    fi
  done
}

if [ "$a_flag" = "true" ]; then
  best_model_metrics $FOLDERS
  echo "Max accuracy of $MAXACC at $MAXACC_SOURCE"
  echo "Max precision of $MAXPREC at $MAXPREC_SOURCE"
  echo "Max recall of $MAXRC at $MAXRC_SOURCE"
else
  get_metrics $FOLDERS
  echo "Acuuracy $ACC, precision $PREC, adjusted recall $ADJ_RECALL"
fi



