#!/usr/bin/env bash

#SBATCH -J DNAMethBS-pipe-methylation-DNA # job name
#SBATCH -q default # partition name
#SBATCH -t 24:00:00 # maximum time of the job in seconds (time limit)
#SBATCH -N 1  # number of parallel tasks
#SBATCH -c 36 # number of cores per task
#SBATCH -o ./slurm_%j_.out # standard output redirection
#SBATCH -e ./slurm_%j.err # standard error redirection



# configure environment with module

module load snakemake/5.11.2
module load rstudio/1.3.1093
module list -t >&2
# ------------------------------------------ Preparing data -----------------------------------------------#

# lanch preparing data script : 

PATHDATA=$(realpath /$1)
CONFIG=$2

# path envirement pipeline
path_work_env='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline'


PATHRESULT=$(cat $path_work_env/config.yml | grep "pathresults" | cut -d "\"" -f2)
mkdir $PATHRESULT/data_symbolic
$path_work_env/script/preparing_symbolic.sh $PATHDATA $PATHRESULT/data_symbolic

# ---------------------------------------------- snakemake -----------------------------------------------#

# launch workflow on node 
# Dry run (simulation)

CORES=36

srun -c $CORES snakemake -n \
--snakefile $path_work_env/snakefile.smk \
--configfile $CONFIG \
--cores $CORES \
--nolock \
--verbose --printshellcmds --reason \

# ---------------------------------------------------------------------------------------------------------#

# Clean all results

$path_work_env/clean.sh



