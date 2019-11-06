train_set=$1

ORDERS=(15 20 25 30)
C_pos=(1 10 100)
TEST=0.3

for order in "${ORDERS[@]}"; do
    for C in "${C_pos[@]}"; do
        qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v ORDER=${order},C=${C},INPUT=${train_set},CPU=10,TEST=${TEST} train-introns.sh
    done
done
