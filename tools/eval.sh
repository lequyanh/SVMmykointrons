 #!/bin/bash

FILES=$1
for f in $FILES/*/
do
  echo "Processing $f file...", $f
  cd $f
  ACC=$(cat *.txt | grep  'Accuracy' | cut -d ' ' -f2)
  TPs=$(cat *.txt | grep  'TP' | cut -d ' ' -f3)
  FPs=$(cat *.txt | grep  'FP' | cut -d ' ' -f3)
  TNs=$(cat *.txt | grep  'TN' | cut -d ' ' -f3)
  FNs=$(cat *.txt | grep  'FN' | cut -d ' ' -f3)

#  PREC=$(( $TPs/($TPs+$FPs) ))
  echo Accuracy $ACC,  $TPs, $FPs, $TNs, $FNs
done
