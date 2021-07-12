
# Pipeline detection DNA methylation : **DNAmethBS** 

  ![bootstrap](Image/logo1.PNG)


A bioinformatics pipeline for the **detection of DNA methylation** from **WGBS** and **MC-Seq** sequencing data. This pipeline, adapted to the needs of the **EpiTree project**, can be fully transposed to other methylation datasets in plants. The pipeline was built as part of the **ANR project EpiTree** ([https://www6.inrae.fr/epitree-project/Le-projet-EPITREE](https://www6.inrae.fr/epitree-project/Le-projet-EPITREE)).

![bootstrap](Image//Le-projet-EPITREE_inra_image.png)

## Authors
* Abdeljalil SENHAJI RACHIK, https://github.com/abdeljalil-senhaji

## Description
This pipeline is based on **Snakemake**. It is able to detect DNA methylation in trees with associated **R** and **python** scripts. 
5 steps (as show in the figure) :

- Cleaning of raw data followed by **quality control**
- **Alignment** with a reference genome
- Elimination of **duplications**
- Detection of methylated cytosines (**mC**) in the three methylation contexts
- Basic **statistical analyzes** on the detection of methylations

`![bootstrap](Images/rg_dag-meth.png)`

## Environment pipeline

![bootstrap](Image/config_pipe.PNG)

#### Pipeline Folder File:
- **Results**: results folder contain :
> - Ordered folder of results for each pipeline job
> - **data_symbolic** : symbolic link of the input data files
>- **Logs**: folder where there is the output of the log files for each sample are used to check if there is an error
>- **Benchmarck**: folder that contains information on resources (snakemake object)
- **Data**: folder where there are the samples (I put 2 samples with 1000bp just for the pipeline test)
- **Rules**: folder that contains each pipeline rule independently (another way to create rules under snakemake) (optionel)
- **script**: folder that contains the scripts created for our pipeline
- **Tools** : tools that are not installed in the cluster
- **DAG**:  visualizing the DAG (Directed Acyclic Graph) of jobs
- **Config_editeur**: scripts for adapting VIM editors to be adapted with snakemake
#### Pipeline Environment File
- ** test_submit_SLURM.sh**: script sbatch of job submission for pipeline test
- **submit_SLURM.sh**: job submission script
- **config.yaml**: configuration file
- **snakefile.smk**: file where there are all the pipeline rules.
- **clean.sh** : script to remove hidden snakemake files and results

submit_SLURM.sh and test_submit_SLURM.sh are specific to the **SLURM cluster**. Therefore, if you are using another cluster, you must create other files.

## System requirements

Here is the list of tool **dependencies** with the version tested and used for the pipeline.

- module load **snakemake** / 5.11.2
- module load **r**/4.0.2
- module load **bioconductor**/3.11 (Installation of package methylKit  v.1.18.0 https://bioconductor.org/packages/release/bioc/html/methylKit.html) 
- module load **pigz** / 2.4 module
- module load **trimgalore** / 0.6.5
- module load **fastqc** / 0.11.9 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- moduel load **multiqc** / 1.9 (https://multiqc.info/)
- module load **samtools**  / 1.11
- module load **bsmapz**

> **Note:** The absolute path of the **BsmapZ** tool (https://github.com/zyndagj/BSMAPz) must be specified in the config.yaml file. if is not installed !!
> 
## Usage

![bootstrap](Image/config_pipe.PNG)

To run the pipeline, it is necessary to fill in these files:

- **Config.yml**: file configuration pipeline 
Modify the necessary pipeline paths : 
```
pathresults: "path/Results/" # Path results pipeline
path-work-env: "path/DNAMethBS_pipeline/" # path envirement pipeline
tools: "path/Tools/BSMAPz-master" #Tools pipeline
REF_GENOME: "/path/.fasta" #reference genome
```
Parameter of the **Trimgalore** tool :
```
# Params_trimgalore
ERROR_RATE: 0.1
LENGTH: 36
QUALITY: 20
ADAPTER1: --illumina
ADAPTER2: AAATCAAAAAAAC
```
> Note : tools are used by **default** except for trimgalore

When using a pipeline on the **MC-Seq** sequencing data, an advantage option in the pipeline allows to specify that the regions **on-target**, for this it is necessary to indicate **1** in the **Target** variable and to put the **file bed** where there are regions of interest in the pipeline **script** folder.
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
Sbatch configuration for the Slurm cluster example:
```
#SBATCH -J DNAMethBS-pipe-methylation-DNA # job name
#SBATCH -q default # partition name
#SBATCH -t 24:00:00 # maximum time of the job in seconds (time limit)
#SBATCH -N 1  # number of parallel tasks
#SBATCH -c 36 # number of cores per task
#SBATCH -o ./slurm_%j.out # standard output redirection
#SBATCH -e ./slurm_%j.err # standard error redirection
```
- **test_submit_SLURM.sh** (Optional): shell script for stimulating pipeline jobs.

> **Note:**
> - It is advisable after each test step to restart the clean.sh script which will delete the hidden **Snakemake** files and the results in order to avoid the problems of submission errors.
> - When launching the pipeline, you have the choice of working on the same pipeline directory or on another location, except it requires each time to change the paths of these files.



#### Step 1: Cloning the git repository
Clone this repository :
```
git clone "https://gitlab................................"
cd DNAmethBS_pipeline
```

#### Step 2: Configure workflow

**>Input :**  

- The pipeline takes **paired end**, R1 and R2 reads as input .
- when launching the pipeline with to submit script (**submit_SLURM.sh**), it starts with a data preparation step with a shell script ([preparing_symbolic.sh  (https://github.com/abdeljalilsenhaji/DNAMethBS_pipeline/blob/main/script/preparing_symbolic.sh "preparing_symbolic.sh")) which allows you to create symbolic links (allows you to assign another access path to a file by pointing to a file name), to adapt the different possible extensions of the reads in **_ R1.fastq.gz \ | _R2. fastq.gz**  (**Note: if there is an error compared to the data, remember to rename the input, indicate the extension for the reads forword and reverse (_R1.fastq.gz \ | _R2.fastq.gz)**.
- The command line for launching the pipeline (Note: **Important to put the entire path to indicate the data and the config file**).

**>Submission command line :**

```
$ sbatch test_submit_SLURM.sh /path/data /path/config.yml
```
```
$ sbatch submit_SLURM.sh /path/data /path/config.yml
```


## Who is maintaining this repository ?

-   Abdeljalil SENHAJI RACHIK- [senhajirachikabdeljalil@gmail.com](senhajirachikabdeljalil@gmail.com) -



















## 1. file configuration
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



## 2. file soumission
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
path_work_env='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline'


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


The user should use the config.yaml file to provide all the necessary inputs for the pipeline:
> the necessary path
- **path-work-env**: working environment
- **pathresults**: environment of the results
- **pathdada**: the data folder path
- **REF_GENOME**: the reference fasta file path
> the parameters and options of the tools:
- ERROR_RATE: 0.1
- LENGTH: 36
- QUALITY: 20
- ADAPT1: --illumina
- ADAPTER2: AAATCAAAAAAAC
> Snakemake options:
- wildcard
- NBSAMPLE:
- THREADS:
> paths to tools:
- tools:


The config.yaml can also be used to modify the parameters (mainly the call of the quality threshold), to define the cluster resources for each program.






