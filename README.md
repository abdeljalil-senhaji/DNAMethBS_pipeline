
# Pipeline detection DNA methylation : **DNAmethBS** 

  ![bootstrap](Image/logo1.PNG)


A bioinformatics pipeline for the **detection of DNA methylation** from **WGBS (Whole Genome Bisulfite Sequencing)** and **MC-Seq (Methyl-capture Sequencing)** sequencing data. This pipeline, adapted to the needs of the **EpiTree project**, can be fully transposed to other methylation datasets in plants. The pipeline was built as part of the **ANR project EpiTree** ([https://www6.inrae.fr/epitree-project/Le-projet-EPITREE](https://www6.inrae.fr/epitree-project/Le-projet-EPITREE)).

![bootstrap](Image//Le-projet-EPITREE_inra_image.png)

## Authors
* Abdeljalil SENHAJI RACHIK, https://github.com/abdeljalil-senhaji

## What is the **DNAmethBS** pipeline?

This pipeline is based on **Snakemake**. It is able to detect DNA methylation for all three methylation contexts (**CG, CHG and CHH**) from **WGBS** and **MC-Seq** sequencing data using bioinformatics tools and scripts **R**, **shell** and **python** associated.
5 steps (as shown in figure):

- Cleaning of raw data followed by **quality control**
- **Alignment** with a reference genome
- Elimination of **duplications**
- Detection of methylated cytosines (**mC**) 
- Extraction methylation DNA in the three methylation contexts (**CG,CHG and CHH**)
- Basic **statistical analyzes** on the detection of methylations

The pipleine is represented in the graph below: 

`![bootstrap](Images/rg_dag-meth.png)`

## Description environment pipeline

```bash

├── clean.sh
├── config.yml
├── DAG
├── Data
├── Editeurs_config
│   ├── configure.sh
│   └── snakemake.vim
├── file_references
│   ├── merge_covered.bed
│   └── Ref_trichocarpa
├── README.md
├── Results
├── Rules
│   ├── bsmapz.smk
│   ├── fastqc.smk
│   ├── methratio.smk
│   ├── multiqc.smk
│   ├── preparing_data_R.smk
│   ├── R_methylkit.smk
│   ├── samtools_markdup.smk
│   ├── samtools_stats.smk
│   ├── splitting_context.smk
│   └── trimgalore.smk
├── script
│   ├── file_renameV2.py
│   ├── preparing_symbolic.sh
│   ├── script_methylkit.R
│   └── splitting.sh
├── snakefile.smk
├── submit_SLURM.sh
├── test_submit_SLURM.sh
└── Tools
    └── BSMAPz-master

```

#### Pipeline environment Folder :
- **Results**: results folder contain :
> - Ordered folder of results for each pipeline job.
> - **data_symbolic** : symbolic link of the input data files.
>- **Logs**: folder where there is the output of the log files for each sample are used to check if there is an error.
>- **Benchmarck**: folder that contains information on resources (snakemake object).
- **Data**: folder where there are the samples (I put 2 samples with 1000bp just for the pipeline test).
- **Rules**: folder that contains each pipeline rule independently (another way to create rules under snakemake) (optionel)
- **script**: folder that contains the scripts created for our pipeline.
- **Tools** : tools that are not installed in the cluster.
- **DAG**:  visualizing the DAG (Directed Acyclic Graph) of jobs.
- **Config_editeur**: scripts for adapting VIM editors to be adapted with snakemake.
- **file_references** : folder where there are the files necessary for the pipeline (for example: the reference fasta genome ...).
#### Pipeline environment File :
- **test_submit_SLURM.sh**: script sbatch of job submission for pipeline test.
- **submit_SLURM.sh**: pipeline submission shell script, This script will be launched as a task on the cluster.
- **config.yaml**: configuration pipeline file.
- **snakefile.smk**: file where there are all the pipeline rules.
- **clean.sh** : script to remove hidden snakemake files and results

submit_SLURM.sh and test_submit_SLURM.sh are specific to the **SLURM cluster**. Therefore, if you are using another cluster, you must create other files.

## System requirements :

Here is the list of tool **dependencies** with the version tested and used for the pipeline.

- module load **snakemake**/5.11.2 (https://snakemake.readthedocs.io/en/stable/)
- module load **r**/4.0.2 (https://cran.r-project.org/bin/windows/base/)
- module load **bioconductor**/3.11 (Installation of package **methylKit** v.1.18.0 https://bioconductor.org/packages/release/bioc/html/methylKit.html) 
- module load **pigz**/2.4 module (https://zlib.net/pigz/)
- module load **trimgalore**/0.6.5 (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- module load **fastqc**/0.11.9 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- moduel load **multiqc**/1.9 (https://multiqc.info/)
- module load **samtools**/1.11 (http://www.htslib.org/)
- module load **bsmapz**/1.1.3 (https://github.com/zyndagj/BSMAPz)

> **Note:** The absolute path of the **BsmapZ** tool (https://github.com/zyndagj/BSMAPz) must be specified in the config.yaml file. if is not installed !!
 
## How to use it ?

Installation and Running the pipeline :

#### Step 1: Cloning the git repository

- Clone the pipeline in your directory :
```
git clone "https://gitlab................................"
cd DNAmethBS_pipeline
```
- Configure config.yml according to the instructions the file.
- Configure submit_SLURM.sh to defene cpu and memory resources for each rule and lunch the workflow.

#### Step 2: Configure workflow

**>Input :**  

- The pipeline takes **paired end**, R1 and R2 reads as input .
- when launching the pipeline with to submit script (**submit_SLURM.sh**), it starts with a data preparation step with a shell script (**preparing_symbolic.sh**) which allows you to create symbolic links (allows you to assign another access path to a file by pointing to a file name), to adapt the different possible extensions of the reads in **_ R1.fastq.gz \ | _R2. fastq.gz**. 
Example of acceptable symbolic links  :    
                                          * (...).1.fastq.gz
                                          * (...).1.fq.gz
                                          * 1_(...).fq.gz
                                          * 1_1_(...).fastq.gz
                                          * _1_1_(...).fq.gz
                                          * (...)_R1.fastq.gz   
                                          * (...).R1.fastq.gz
                                          * Ps. (...) Corresponds to the baseName of sample
**Note: whatever is written in upper or lower case is acceptable !!                                     
**Note: if there is an error compared to the data, remember to rename the input, indicate the extension for the reads forword and reverse (_R1.fastq.gz \ | _R2.fastq.gz)**.
- The command line for launching the pipeline (Note: **Important to put the entire path to indicate the data and the config file**).

## Usage

![bootstrap](Image/config_pipe.PNG)

To run the pipeline, it is necessary to fill in these files:

- **Config.yml**: file configuration pipeline 
Modify the necessary pipeline paths : 
```
pathresults: "path/Results/"                                   # Path results pipeline
path-work-env: "path/DNAMethBS_pipeline/"                      # path envirement pipeline
tools: "path/Tools/BSMAPz-master"                              # Tools pipeline
REF_GENOME: "/path/.fasta"                                     # reference genome
```

Parameter of the **Trimgalore** tool :
```
# Params_trimgalore
ERROR_RATE: 0.1                             # Maximum allowed error rate (no. of errors divided by the length of the matching region)
LENGTH: 36                                  # Discard reads that became shorter than length INT because of either quality or adapter trimming. 
QUALITY: 20                                 # Trim low-quality ends from reads in addition to adapter removal.
ADAPTER1: --illumina                        # Adapter sequence to be trimmed. If not specified explicitly
ADAPTER2: AAATCAAAAAAAC                     # Optional adapter sequence to be trimmed off read 2 of paired-end files.
```
> Note : tools are used by **default** except for trimgalore
When using a pipeline on the sequencing data **MC-Seq**, an advantage in the pipeline allows to specify that the regions **On-target**, for this it is necessary to indicate **1** in the **Target** and place the **file bed** where there are regions of interest in the pipeline's **script** folder. If you do not use the MC-seq technique, indicate **0**
```
# SELECT EITHER BY SPECIFY THE ONTARGET (1) OR NOT (0)
Target: 0
#Path the file bed 
pathbed: "/path/script/merge_covered.bed"
```
- **submit_SLURM.sh**: pipeline submission shell script, This script will be launched as a task on the cluster, acting as a "master" script executing all the tasks required at all stages of the pipeline.
```
path_work_env= '/path/envirennement_pipeline'
```
Sbatch configuration for the **Slurm cluster** example:
```
#SBATCH -J DNAMethBS-pipe-methylation-DNA                       # job name
#SBATCH -q default                                              # partition name
#SBATCH -t 24:00:00                                             # maximum time of the job in seconds (time limit)
#SBATCH -N 1                                                    # number of parallel tasks
#SBATCH -c 36                                                   # number of cores per task
#SBATCH -o ./slurm_%j.out                                       # standard output redirection
#SBATCH -e ./slurm_%j.err                                       # standard error redirection
```
- **test_submit_SLURM.sh** (Optional): shell script for stimulating pipeline jobs.

> **Note:**
> - It is advisable after each test step to restart the clean.sh script which will delete the hidden **Snakemake** files and the results in order to avoid the problems of submission errors.
> - When launching the pipeline, you have the choice of working on the same pipeline directory or on another location, except it requires each time to change the paths of these files.

#### Launch the workflow

**>Submission command line :**

You can start by boing a dry run test to visualize all the jobs and commands:
```
$ sbatch test_submit_SLURM.sh /path/data/ /path/config.yml
```
To run on slurm:
```
$ sbatch submit_SLURM.sh /path/data/ /path/config.yml
```

## 1. Exemple file configuration
```

######################### Configuration snakemake #############################

## -----------------------------------> variable path pipeline 

# Path of results

pathresults: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/Results/"

# Switch reference genome Populus Trichocarpa V3/V4 

REF_GENOME: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/Ref_trichocarpa/Ptrichocarpa_444_v3.0.fa"
#REF_GENOME: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/Ref_trichocarpa/Ptrichocarpa_533_v4.0.fa"

# Number of samples for modeling THREADS of each rule in the pipeline

NBSAMPLE: 1


## -----------------------------------> Settings softwars


# All default settings except trimgalore :

# Params_trimgalore
 
ERROR_RATE: 0.1
LENGTH: 36
QUALITY: 20
ADAPTER1: --illumina
ADAPTER2: AAATCAAAAAAAC

## rule duplicate specify ON-Target


# SELECT EITHER BY SPECIFY THE ONTARGET (1) OR NOT (0)
Target: 0

#Path the file bed 
pathbed: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/script/merge_covered.bed"

## -----------------------------------> path pipeline 

#  Path environnement pipeline 

path-work-env: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/"

# Software_path GitHub Repository (https://github.com/zyndagj/BSMAPz)

tools: "/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/Tools/BSMAPz-master"

# Wildcard 

CONTEXT: ['CG', 'CHH', 'CHG']
READS: ['1', '2']

#--------------------------------------------------------------------------------------------
```



## 2. Exemple file soumission
```
#!/usr/bin/env bash

#SBATCH -J DNAMethBS-pipe-methylation-DNA # job name
#SBATCH -q default # partition name
#SBATCH -t 24:00:00 # maximum time of the job in seconds (time limit)
#SBATCH -N 1  # number of parallel tasks
#SBATCH -c 36 # number of cores per task
#SBATCH -o ./slurm_%j.out # standard output redirection
#SBATCH -e ./slurm_%j.err # standard error redirection

# Dry run (simulation)

# configure environment with module

module load snakemake/5.11.2
module load rstudio/1.3.1093
module list -t >&2
# ------------------------------------------ Preparing data -----------------------------------------------#

# lanch preparing data script : 

PATHDATA=$(realpath $1)
CONFIG=$2

# path envirement pipeline
path_work_env='/path/DNAMethBS_pipeline'


PATHRESULT=$(cat $CONFIG | grep "pathresults" | cut -d "\"" -f2)
mkdir $PATHRESULT/data_symbolic
$path_work_env/script/preparing_symbolic.sh $PATHDATA $PATHRESULT/data_symbolic

# ---------------------------------------------- snakemake -----------------------------------------------#

# launch workflow on node
CORES=36

srun -c $CORES snakemake \
--snakefile $path_work_env/snakefile.smk \
--configfile $CONFIG \
--cores $CORES \
--nolock \
--verbose --printshellcmds --reason \

#--forceall
# ---------------------------------------------------------------------------------------------------------#
```





## Who is maintaining this repository ?

-   Abdeljalil SENHAJI RACHIK -[senhajirachikabdeljalil@gmail.com](senhajirachikabdeljalil@gmail.com)-




