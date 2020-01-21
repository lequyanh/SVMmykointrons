INTRON_TRAIN_FILE="intron-train-dataset.csv"
ROOT="lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh"
NCPUS=12

division=ascomycota

echo "Assuming, the file ${INTRON_TRAIN_FILE} exist at ${ROOT}/data/intron, run this command on metacentrum"
echo "bash grid_search_intron.sh $INTRON_TRAIN_FILE"

# ------------------------------------- wait for gridsearch to finish... -----------------------------------------

# this will produce model files, which will be downloaded
models_folder="${division}/intron/models"
mkdir -p ./$models_folder
scp ${ROOT}/results/train-introns-* $models_folder

# and evaluate their performance on validation dataset. Note the parameters of the best performing one
bash validate_grid_models_introns.sh "../preprocessing/intron/${division}/intron_candidates/intron-test-dataset.csv" $models_folder "grid_validation_results/${division}/intron"

# suppose those are the best parameters
C_pos=200
C_neg=1
order=6

# Train the overall intron model with them using all data
$PYTHON ../classification/train-introns.py "../preprocessing/intron/${division}/intron_candidates/$INTRON_TRAIN_FILE" ${order} ${C_neg} ${C_pos} -c ${NCPUS}
