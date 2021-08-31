# Enable conda activate from shell
eval "$(conda shell.bash hook)"

# INSTALL TENSORFLOW-KERAS ENVIRONMENT
conda create -n mykointron python=3.7.4 --yes
conda activate mykointron

conda install pandas --yes
conda install -c anaconda biopython --yes
conda install docopt --yes
conda install scikit-learn --yes
conda install -c anaconda keras --yes

pip3 install --upgrade tensorflow

# INSTALL SHOGUN ENVIRONMENT
conda create -n mykointron_shogun python=3.6 --yes
conda activate mykointron_shogun

conda install pandas --yes
conda install -c anaconda biopython --yes
conda install docopt --yes
conda install scikit-learn --yes
conda install -c conda-forge shogun --yes
