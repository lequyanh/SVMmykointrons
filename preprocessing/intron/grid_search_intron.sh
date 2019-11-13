train_set=$1

ORDERS=(3 4 5)
C_pos=(1 2 3 4)
TEST=0.2

for order in "${ORDERS[@]}"; do
    for C in "${C_pos[@]}"; do
        # qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v ORDER=${order},C=${C},INPUT=${train_set},CPU=10,TEST=${TEST} train-introns.sh
        python ../../classification/train-introns.py "$train_set" "$order" 1 "$C" -t $TEST -c 10
    done
done
