Overview
==========

![Diagram](/diagram.png)

We present a pipeline for removing fungal introns from meta-genome assemblies. The output is the original assembly
cleansed of potential intron sequences. 

The pipeline is separated into 3 stages, each represented by a bash script.

Breaking down the steps into separate stages is necessary due to asynchronous and distributed nature of assembly processing.
The stages are following:
1) *Split the assembly into smaller batches*
   * The size of the batches should be less than 70Mb for reasonable processing times
2) *Process the batches (in parallel)*
   * Here typically the batches are submitted as a job for execution on some cluster
   * Small assemblies can be processed locally
3) *Assemble the results*
    * Once all partial results are complete, merge them into a single cleaned assembly
    * A CSV with removed intron sequences is also included

Installation
============
1. Run the `install.sh` script. 
   
   a. If your shell does not permit running conda from a script, execute the commands manually:
```
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
├── assembly_batches/
│   ├── ctg_k141_2751450_0.fa
│   ├── ctg_k141_2751450_1.fa
│   └── ...
├── ctg_k141_2751450_no_duplicates.fa
└── ctg_k141_2751450.fa
```
The `assembly_batches` directory simply contains the batches themselves. The file suffix signals the index of the first sequence in the batch
Note that the sharding process implicitly removes duplicated records from the assembly. 
The pipeline will later work only with the non-duplicated sequences.

2. After the second step (**independent intron removal from each batch**). For simplicity, we only chose to process the
negative strand (the `minus` suffix of files)
```
.
├── assembly_batches/
├── results/
│   ├── ctg_k141_2751450_0.fa_results_minus/
│   ├── ctg_k141_2751450_1.fa_results_minus/
├── scratchdir/
├── assembly_batches.txt
├── ctg_k141_2751450_no_duplicates.fa
└── ctg_k141_2751450.fa
```
Here a few new items appear. The most important is the `results` directory with partial results for each batch.
`scratchdir` contains only temporary processing files and is of no importance apart from debugging. 
`assembly_batches.txt` holds the list of batch names to process. In the basic scenario, all batches are present; 
In specific cases, where the user needs only a subset of batches processed (e.g. when some cluster node crashes), 
he or she can modify this file to specify the subset required.

3. After the final step (**combining partial results**).
```
.
├── assembly_batches/
├── results/
├── scratchdir/
├── assembly_batches.txt
├── ctg_k141_2751450_no_duplicates.fa
├── ctg_k141_2751450.fa
├── cut-coords-minus-strand-full.csv
└── pruned_ctg_k141_2751450_minus.fa
```
Two new files appear - the pruned assembly (with non-duplicated sequences), and the CSV with cut coordinates.
Both files are assembled from the partial results inside the `results` folder. End of processing

### IMPORTANT: Distributed job submission
Though the entire pipeline can be executed locally for small assemblies, the main use-case is to send individual batches 
as distributed tasks to some high-performance cluster. We are aware that each cluster has a different infrastructure
and different job submission policy. This particular part of the pipeline must be therefore left for the user to fill in.

For convenience, we isolated the submission section into the `submit_job()` functions inside the `batch_pipeline.sh` script.
Basically, the user only needs to modify a single line that specifies which script and with which arguments should be executed.
The script will have had already its working directory with all necessary files prepared. 

Quickstart
==========
Within the `./pipeline` folder, run:

1) **Split the assembly into smaller batches for parallel/distributed processing:**

`bash shard_assembly.sh -p project_path -n seqs_per_batch`

2) **Process the batches:**

`bash batch_pipeline.sh -m models_settings -p project_path -s strand`

3) **Combine the results**

`bash combine_results.sh -p project_path`

### Example 1 (demo)
As a demonstration, we will process a test contig ctg_k141_2751450. The contig is part of the sources and is located
in the test directory *test/projects/project_ctg_k141_2751450*.

First, split the assembly into batches with 1 sequence each. The batches will be saved into the project directory
under the folder `assembly_batches`

`bash shard_assembly.sh -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450/ -n 1`

Next, we process the batches, removing introns from both strands. We will use SVM models trained on Basidiomycota species as indicated by 
the `svmb` argument. The computations will be performed locally:

`bash batch_pipeline.sh -m svmb -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450 -s both`

Finally, we combine the partial results into a single assembly purged from introns. As we were processing both strands, 
the call results in two separate fastas one for each strand: `pruned_ctg_k141_2751450_no_duplicates_minus.fa` and `pruned_ctg_k141_2751450_no_duplicates_plus.fa`

`bash combine_results.sh -p /home/john/mycointrons/test/projects/project_ctg_k141_2751450`

### Example 2 (*Kocim1*; scaffold_1)

The test directory also contains a more practical example.  
Now, we will be removing introns from the first scaffold of the *Kocim1* fungi. The process is identical to the one above:

```
bash shard_assembly.sh -p /home/john/mycointrons/test/projects/project_Kocim1_reduced/ -n 1

bash batch_pipeline.sh -m svmb -p /home/john/mycointrons/test/projects/project_Kocim1_reduced -s both

bash combine_results.sh -p /home/john/mycointrons/test/projects/project_Kocim1_reduced
```

The process should take approx. 1h for each strand. 

We also provide reference results stored in the `ref/` subdirectory of the project for verification. 

Details
=====
For model setting, use one of the following
* 'svmb' (standing for SVM models trained on Basidiomycota)
* 'random' (dummy classification; placeholder for future models)

The process of intron removal is roughly:
1) Find all donor dimers and perform splice site classification
2) Find all acceptor dimers, remove orphan candidates (AG with no GT in acceptable range) and perform splice site classification
3) Pair positively classified splice site candidates to form an intron candidate dataset
4) Classify the intron dataset
5) Cut positively classified introns. Overlaps are resolved with length prior distribution cut-off
