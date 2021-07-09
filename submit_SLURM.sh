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


PATHRESULT=$(cat CONFIG | grep "pathresults" | cut -d "\"" -f2)
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

# Useful information to print

#echo '########################################'
#echo 'Date:' $(date --iso-8601=seconds)
#echo 'User:' $USER
#echo 'Host:' $HOSTNAME
#echo 'Job Name:' $SLURM_JOB_NAME
#echo 'Job ID:' $SLURM_JOB_ID
#echo 'Array task ID:' ${SLURM_ARRAY_TASK_ID}
#echo 'Number of nodes assigned to job:' $SLURM_JOB_NUM_NODES
#echo 'Total number of cores for job (?):' $SLURM_NTASKS
#echo 'Number of requested cores per node:' $SLURM_NTASKS_PER_NODE
#echo 'Nodes assigned to job:' $SLURM_JOB_NODELIST
#echo 'Directory:' $(pwd)
## Detail Information:
#echo 'scontrol show job:'
#scontrol show job $SLURM_JOB_ID
#echo '########################################'

