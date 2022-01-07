# Enable conda activate from shell
eval "$(conda shell.bash hook)"

conda create -n mykointron_shogun python=3.6 --yes
conda activate mykointron_shogun

conda install pandas --yes
conda install -c anaconda biopython --yes
conda install docopt --yes
conda install scikit-learn --yes
conda install -c conda-forge shogun --yes

# Return to the base environment
conda activate