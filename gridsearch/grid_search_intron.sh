train_set=$1

ORDERS=(4 5 6 7 8)
C_pos=(3 4 5 6)
TEST=0.2

for order in "${ORDERS[@]}"; do
    for C in "${C_pos[@]}"; do
        # qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v ORDER=${order},C_POS=${C},INPUT=${train_set},CPU=10,TEST=${TEST} train-introns.sh
        python ../classification/train-introns.py "${train_set}" "${order}" 1 "${C}" -t ${TEST} -c 12 -o "intron-model-d-${order}-Cpos-${C}.hd5"

        resultdir="intron-model-${train_set}-d-${order}-Cpos-${C}"
        mkdir "${resultdir}"

        mv train-introns*.log "${resultdir}"
        mv intron-model*.hd5 "${resultdir}"
    done
done
