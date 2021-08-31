conda create -n mykointron python=3.7.4 --yes
conda activate mykointron

conda install pandas --yes
conda install -c anaconda biopython --yes
conda install docopt --yes
conda install scikit-learn --yes
conda install -c anaconda keras --yes

pip3 install --upgrade tensorflow