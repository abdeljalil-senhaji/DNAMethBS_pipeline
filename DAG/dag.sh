#!/usr/bin/env bash

# Generate the dag files

module load snakemake/5.11.2

WORKFLOW='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/snakefile.smk'
CONFIG='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/config.yml'
OUTPATH='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/DAG'

# With samples

snakemake \
--snakefile $WORKFLOW \
--configfile $CONFIG \
--dag > $OUTPATH/dag.dot
dot -Tpng $OUTPATH/dag.dot > $OUTPATH/dag-meth.png

# Only rules

snakemake \
--snakefile $WORKFLOW \
--configfile $CONFIG \
--rulegraph > $OUTPATH/rg_dag.dot
dot -Tpng $OUTPATH/rg_dag.dot > $OUTPATH/rg_dag-meth.png

# remove intermediate files
rm $OUTPATH/dag.dot $OUTPATH/rg_dag.dot


