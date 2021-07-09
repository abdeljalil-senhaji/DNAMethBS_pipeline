#!/usr/bin/env bash

## clean

module load snakemake/5.11.2

RESULTSPATH='/env/cng/proj/LEE/pipelines/DNAMethBS_pipeline/Results'

core=1

snakemake -s snakefile.smk --cores $core --delete-all-output --dry-run
snakemake -s snakefile.smk --cores $core --delete-all-output

rm -r $RESULTSPATH/*
rm -r .snakemake/
