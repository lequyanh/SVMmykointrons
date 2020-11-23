USAGE
=====

The call is:
`bash main.sh PATH_TO_METAGENOM MODEL_SETTINGS`

for model setting, use one of the following
    * 'svmb' (standing for SVM models trained on Basidiomycota; very slow)
    * 'nn100' (standing for neural net with 0-100 windows; faster)
    * 'nn200' (standing for neural net with 200-200 windows; slower, more accurate)

The process is roughly:
    1) Split metagenom file into smaller chunks. This allows more parallel and distributed approach
    2) Perform intron cutting
        a) Find all donor dimers and perform splice site classification
        b) Find all acceptor dimers, remove orphan candidates (AG with no GT in acceptable range) and perform splice site classification
        c) Pair positively classified splice site candidates to form an intron candidate dataset
        d) Classify the intron dataset
        e) Cut positively classified introns. Overlaps are resolved with length prior distribution cut-off
    3) Cuts annotation

The process can be tweeked - see below for details on how to perform partial steps of the pipeline

CONTENTS
========

(0) main.sh
(1) annotation/
(2) classification/
    (a) classify-introns.py
    (b) classify-splice-sites.py
        python classify-splice-sites.py ../data/train/donor/Exova1-donor-windows.csv ../gridsearch/bestmodels/donor_model.hd5 70 70 donor -v -c 10
    (c) train-introns.py
    (d) train-splice-sites.py
(3) pipeline/
    (a) pipeline.sh
    (b) batch_pipeline.sh

(4) taxonomy.csv
(5) tools/
    (a) extend-fasta.py
    (b) fastalib.py
    (c) generate-pairs.py
    (d) process-gff.py


(0) main.sh

(1) annotation

(2) classification

    Contains scripts for model training and classification (for both introns & splice sites).
    The scripts depend on Shogun and Pandas (see section INSTALLATION).

    The scripts provide command line interfaces. Use '-h' option to display help.

    Inputs of the train scripts are expected to be CSV files with columns 'sequence' and 'label'.
    Inputs of the classify scripts are expected to be CSV files with column 'sequence'.


(3) pipeline

    Contains an implementation of the classification process.
    The main component is a bash script called `pipeline` and 'batch_pipeline' respectively
    The script calls other python scripts and performs the whole process.

    This is the main part of intron classification process. Can be used on its own but requires many parameters.
    On the other hand can allow processing of a single assembly and also provides cuts validation in case true introns are known


(4) taxonomy.csv

    Contains table with taxonomical classification.


(5) tools

    extend-fasta.py extends the input fasta with scaffolds representing negative strands
    of the original scaffolds

    fastalib.py is a small library defining convenient functions for processing of FASTA files

    generate-pairs.py generates intron candidates by pairing donors and acceptors

    process-gff.py parses a GFF file and finds introns
    TIP: output of this script has the same syntax as the input of extract-fasta
         it is therefore suitable for piping

    See comments inside the scripts for more information.



INSTALLATION
============

(1) PYTHON libraries
      All the scripts assume Python3 is used.

      We need two separate environments with similar libraries, since Shogun and Keras don't work together in one env.

        conda create -n mykointron_shogun python=3.6
        conda config --env --add channels conda-forge anaconda
        conda activate mykointron_shogun
        conda install pandas
        conda install -c anaconda biopython
        conda install docopt
        conda install scikit-learn
        (pip install gffutils)
        conda install -c conda-forge shogun

        conda create -n mykointron python=3.6
        conda activate mykointron
        conda install pandas
        conda install -c anaconda biopython
        conda install docopt
        conda install scikit-learn
        conda install -c anaconda keras



