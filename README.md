Installation
============
1. Run the `install.sh` script. 
   
   a. If your shell does not permit running conda from a script, execute the commands manually:
```
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

# Return to the base environment
conda activate
```

2. Download the models from TODO and place them to `pipeline/models/`

Overview
==========
We present a pipeline for removing fungal introns from meta-genome assemblies. The output is the original assembly
purged of potential intron sequences. 

The pipeline is separated into 3 stages, each represented by a bash script.

Breaking down the steps into separate stages is necessary due to asynchronous and distributed nature of assembly processing.
The stages are following:
1) *Split the assembly into smaller shards*
   * The size of the shards should be less than 70Mb for reasonable processing times
2) *Process the shards (in parallel)*
   * Here typically the shards are submitted as a job for execution on some cluster
   * Small assemblies can be processed locally
3) *Combine the results*
    * Once the partial results are complete, merge them into a single cleaned assembly

### Project folders    
The asynchronous and distributed nature of the task requires a central location where partial results are stored.
We therefore introduce a concept of a *project folder*.

At the beginning, the folder only contains the assembly fasta to process. As the pipeline progresses, the folder gets
populated by partial results and helper files. 
We will now take a look at the structure of the folder at each stage of the pipeline (using a test assembly *ctg_k141_2751450.fa* ). 

The starting appearance is simply:

```
.
└── ctg_k141_2751450.fa
```

1. After the first step (**assembly sharding**). :
```
.
├── assembly_shards/
│   ├── ctg_k141_2751450_0.fa
│   ├── ctg_k141_2751450_1.fa
│   └── ...
├── ctg_k141_2751450_no_duplicates.fa
└── ctg_k141_2751450.fa
```
The `assembly_shards` directory simply contains the shards themselves. The file suffix signals the index of the sequence first present in the shard
Note that the sharding process implicitly removes duplicated records from the assembly. 
The pipeline will later work only with the non-duplicated sequences.

2. After the second step (**independent intron removal from each shard**). For simplicity, we only chose to process the
negative strand (the `minus` suffix of files)
```
.
├── assembly_shards/
├── results/
│   ├── ctg_k141_2751450_0.fa_results_minus/
│   ├── ctg_k141_2751450_1.fa_results_minus/
├── scratchdir/
├── assembly_shards.txt
├── ctg_k141_2751450_no_duplicates.fa
└── ctg_k141_2751450.fa
```
Here a few new items appear. The most important is the `results` directory with partial results for each shard.
`scratchdir` contains only temporary processing files and is of no importance apart from debugging. 
`assebmly_shards.txt` holds the list of shard names to process. In the basic scenario all shards are present; 
In specific cases, where the user needs only a subset of shards processed (e.g. when some cluster node crashes), 
he or she can modify this file to specify the subset required.

3. After the final step (**combining partial results**).
```
.
├── assembly_shards/
├── results/
├── scratchdir/
├── assembly_shards.txt
├── ctg_k141_2751450_no_duplicates.fa
├── ctg_k141_2751450.fa
├── cut-coords-minus-strand-full.csv
└── pruned_ctg_k141_2751450_minus.fa
```
Two new files appear - the pruned assembly (with non-duplicated sequences), and the CSV with cut coordinates.
Both files are assembled from the partial results inside the `results` folder. End of processing

### IMPORTANT: Distributed job submission
Though the entire pipeline can be executed locally for small assemblies, the main use-case is to send individual shards 
as a distributed task to some high-performance cluster. We are aware that each cluster has a different infrastructure
and different job submission policy. This particular part of the pipeline must be therefore left for the user to fill in.

For convenience, we isolated the submission section into the `submit_job()` functions inside the `batch_pipeline.sh` script.
Basically, the user only needs to modify a single line that specifies which script and with which arguments should be executed.
The script will have had already its working directory with all necessary files prepared. 

Quickstart
==========
Within the `./pipeline` folder, run:

1) **Split the assembly into smaller shards for parallel/distributed processing:**

`bash shard_assembly.sh -p project_path -n seqs_per_shard`

2) **Process the shards:**

`bash batch_pipeline.sh -m models_settings -p project_path -s strand`

3) **Combine the results**

`bash combine_results.sh -p project_path`


As a demonstration, we will process a test contig ctg_k141_2751450. The contig is part of the sources and is located
in the test directory *test/projects/project_ctg_k141_2751450*.

First, split the assembly into shards with 1 sequence each. The shards will be saved into the project directory
under the folder `assembly_shards`

`bash shard_assembly.sh -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450/ -n 1`

Next, we process the shards, removing introns from both strands. We will use SVM models trained on Basidiomycota species as indicated by 
the `svmb` argument. The computations will be performed locally:

`bash batch_pipeline.sh -m svmb -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450 -s both`

Finally, we combine the partial results into a single assembly purged from introns. As we were processing both strands, 
the call results in two separate fastas one for each strand (named `pruned_ctg_k141_2751450_no_duplicates_minus.fa` and `pruned_ctg_k141_2751450_no_duplicates_plus.fa`)

`bash combine_results.sh -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450`

Details
=====

for model setting, use one of the following
* 'svmb' (standing for SVM models trained on Basidiomycota; very slow)
* 'random' (dummy classification; placeholder for future models)

The process of intron removal is roughly:
1) Find all donor dimers and perform splice site classification
2) Find all acceptor dimers, remove orphan candidates (AG with no GT in acceptable range) and perform splice site classification
3) Pair positively classified splice site candidates to form an intron candidate dataset
4) Classify the intron dataset (only for SVM models)
5) Cut positively classified introns. Overlaps are resolved with length prior distribution cut-off
