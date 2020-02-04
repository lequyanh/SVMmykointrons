train_set=$1

ORDERS=(4 5 6 7 8)
C_pos=(3 4 5 6 7)
TEST=0.2

for order in "${ORDERS[@]}"; do
    for C in "${C_pos[@]}"; do
        qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v ORDER="${order}",C_POS="${C}",INPUT="${train_set}",CPU=10,TEST=${TEST} train_introns_meta.sh
    done
done
