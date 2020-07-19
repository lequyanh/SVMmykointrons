CONTENTS
========

(1) extract-fasta/
(2) classification/
    (a) classify-introns.py
    (b) classify-splice-sites.py
        python classify-splice-sites.py ../data/train/donor/Exova1-donor-windows.csv ../gridsearch/bestmodels/donor_model.hd5 70 70 donor -v -c 10
    (c) train-introns.py
    (d) train-splice-sites.py
(3) pipeline/
        bash pipeline ~/Desktop/mykointrons-data/data/Assembly/Exova1_AssemblyScaffolds.fasta ../gridsearch/bestmodels/donor_model.hd5 ../gridsearch/bestmodels/acceptor_model.hd5 intron-model-Cpos-3-deg-5.hd5

(4) taxonomy.csv
(5) tools/
    (a) extend-fasta.py
    (b) fastalib.py
    (c) generate-pairs.py
    (d) process-gff.py



(1) extract-fasta

    Tool that extracts subsequences from FASTA file.

    Use `make` to build a binary.
    Uses `g++` for the build (tested with ver. 7.4.0)

    Use `./extract-fasta -h` to understand the usage.


(2) classification

    Contains scripts for model training and classification (for both introns & splice sites).
    The scripts depend on Shogun and Pandas (see section INSTALLATION).

    The scripts provide command line interfaces. Use '-h' option to display help.

    Inputs of the train scripts are expected to be CSV files with columns 'sequence' and 'label'.
    Inputs of the classify scripts are expected to be CSV files with column 'sequence'.


(3) pipeline

    Contains a prototype implementation of the classification process.
    The main component is a bash script called `pipeline`.
    The script calls other python scripts and performs the whole process.
    The script contains detailed comments.


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

  All the scripts assume Python3 is used.

  Libraries Shogun and Pandas are required.
  Instructions on how to install
    - Shogun are at http://shogun.ml/install
    - Pandas are at https://pandas.pydata.org/pandas-docs/stable/install.html

    conda create -n mykointron python=3.6
    conda config --env --add channels conda-forge anaconda
    conda install pandas
    conda install -c conda-forge shogun
    conda install -c anaconda keras
    conda install -c anaconda biopython
    (pip install gffutils)
