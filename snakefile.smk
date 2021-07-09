# coding: utf-8
###################################################################################################
#                                               DNAMethBS                     		          #
#				     PIPELINE METHYLATION BISULFITE DNA 		          #
###################################################################################################


"""
__Author__ = "SENHAJI RACHIK Abdeljalil"
__Version__ = "1.0.0"
__Description__: Bioinformatics pipeline for detecting DNA methylation.
__Project__: ANR EPITREE
__Date__: LAST UPDATE >> 30-06-2021

"""


#---------------------------------------pipeline settings ------------------------------------------#

#-------------------> Before making a modification, be sure to read the Readme <--------------------#

#----------------------------------Library needed for snakemake ------------------------------------#

import os
from os.path import join
import glob
import yaml
from snakemake.exceptions import MissingInputException
import collections
import shutil
from pathlib import Path
import re
import stat
import time
import datetime
import json
import copy
import functools
import subprocess as sp
from itertools import product, chain
from contextlib import contextmanager
import string
import asyncio

# ------------------------------------- Setup Snakefile ---------------------------------------------#

shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")
localrules : all

#Path to the config file, defined parameters like path of first input files
#configfile: "./config.yml"

# --------------------Choice of the condition with extraction of ontarget or not --------------------#

if config['Target'] == 1:
    ruleorder: samtools_markdup_ontarget > samtools_markdup
else:
    ruleorder: samtools_markdup > samtools_markdup_ontarget


# --------------------- Déclaring some variable required for the subsequent steps -------------------#

#To get the list of sample files to process I used the glob_wildcards() function.

pathdada = config['pathresults'] + 'data_symbolic/'

SAMPLES, = glob_wildcards(pathdada + '{sample}_R1.fastq.gz')

#############################################################################################
#                                    pipeline rules                                         #
#############################################################################################

rule all:
 input:
  expand(config['pathresults'] + '01_trimgalore/{sample}/{sample}_R1_val_1.fq.gz', sample=SAMPLES),
  expand(config['pathresults'] + '01_trimgalore/{sample}/{sample}_R2_val_2.fq.gz', sample=SAMPLES),      
  expand(config['pathresults'] + '01_trimgalore/{sample}/{sample}_R1.fastq.gz_trimming_report.txt', sample=SAMPLES),
  expand(config['pathresults'] + '01_trimgalore/{sample}/{sample}_R2.fastq.gz_trimming_report.txt', sample=SAMPLES),
  expand(config['pathresults'] + '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.html', sample=SAMPLES, reads=config['READS']),
  expand(config['pathresults'] + '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.zip', sample=SAMPLES, reads=config['READS']),
  config['pathresults'] + '03_multiqc',
  expand(config['pathresults'] + '04_bsmapz/{sample}/mapped_{sample}.bam', sample=SAMPLES),
  expand(config['pathresults'] + '05_samtools_stats/{sample}/statistics_{sample}.txt', sample=SAMPLES),
  expand(config['pathresults'] + '06_markdup/{sample}/{sample}.dedup.bam', sample=SAMPLES),
  expand(config['pathresults'] + '07_methratio/{sample}/meth_{sample}.txt', sample=SAMPLES),
  expand(config['pathresults'] + '08_splitting_context/{sample}/{context}-{sample}.txt', sample=SAMPLES, context=config['CONTEXT']),
  expand(config['pathresults'] + '09_preparing_data_R/{context}/{context}-{sample}.txt', sample=SAMPLES, context=config['CONTEXT']),
  expand(config['pathresults'] + '10_methylkit/{context}/{context}.csv', context=config['CONTEXT']),
  expand(config['pathresults'] + '10_methylkit/{context}/{context}.pdf', context=config['CONTEXT'])

############################ Rule trimgalore ################################################
# 		Step 1 : Cleaning of raw data followed by quality control		    #
############################ Rule trimgalore ################################################

rule trimgalore:
 input:
  fwd = pathdada + '{sample}_R1.fastq.gz',
  rev = pathdada + '{sample}_R2.fastq.gz'
 output:
  out1 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R1_val_1.fq.gz',
  out2 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R2_val_2.fq.gz',
  temp1 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R1.fastq.gz_trimming_report.txt',
  temp2 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R2.fastq.gz_trimming_report.txt'
 message:
   "Trimming reads {input.fwd} and {input.rev}"
 params:
  fq = os.path.join(config['pathresults'] + '01_trimgalore/{sample}/')
 threads:
  workflow.cores/(config['NBSAMPLE']*1)
 log:
  out = config['pathresults'] + 'logs/01_trimgalore/{sample}/trimmed_{sample}.o',
  err = config['pathresults'] + 'logs/01_trimgalore/{sample}/trimmed_{sample}.e',
  rc =  config['pathresults'] + 'logs/01_trimgalore/{sample}/trimmid_{sample}.rc'
 benchmark:
   config['pathresults'] + 'benchmarks/01_trimgalore/{sample}/trimmed_{sample}.tsv'
 shell:
  '''
   module load trimgalore/0.6.5
   module load pigz/2.4
   (
    trim_galore {input.fwd} {input.rev} --paired {config[ADAPTER1]} -a2 {config[ADAPTER2]} -o {params.fq} --gzip -j {threads} 
   ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

############################################## FASTQC #####################################
# 				Step 2 : quality control reads				  #
###########################################################################################

rule fastqc:
 input:
  inp = os.path.join(config['pathresults'], '01_trimgalore/{sample}/{sample}_R{reads}_val_{reads}.fq.gz')
 output:
  os.path.join(config['pathresults'], '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.html'),
  os.path.join(config['pathresults'], '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.zip')
 threads: 
  workflow.cores/(config['NBSAMPLE']*2)
 message:
  'Quality check of raw data with fastqc.'
 params:
  reg = os.path.join(config['pathresults'] + '02_fastqc/{sample}/')
 log:
  out = config['pathresults'] + 'logs/02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.o',
  err = config['pathresults'] + 'logs/02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.e',
  rc =  config['pathresults'] + 'logs/02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.rc'
 benchmark:
   config['pathresults'] + 'benchmarks/02_fastqc/{sample}_R{reads}_val_{reads}_fastqc.tsv'
 shell:
  '''
  module load fastqc/0.11.9
  (
   fastqc {input.inp} -o {params.reg} -t {threads}
  ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

########################################## multiQC  #########################################
# 			step 3 : to gather all results in an html file			    #
#############################################################################################

rule multiqc:
 input:
  expand(config['pathresults'] + '01_trimgalore/{sample}/{sample}_R{reads}.fastq.gz_trimming_report.txt', sample=SAMPLES, reads=config['READS']),
  expand(config['pathresults'] + '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.html', sample=SAMPLES, reads=config['READS']),
  expand(config['pathresults'] + '02_fastqc/{sample}/{sample}_R{reads}_val_{reads}_fastqc.zip', sample=SAMPLES, reads=config['READS']),
  expand(config['pathresults'] + '05_samtools_stats/{sample}/statistics_{sample}.txt', sample=SAMPLES)
 output:
  directory(config['pathresults'] + '03_multiqc')
 message:
  'Quality check of raw data with multiqc.'
 threads: workflow.cores/(config['NBSAMPLE']*1)
 benchmark:
   config['pathresults'] + 'benchmarks/03_multiqc/multiqc_report.tsv'
 log:
  out = config['pathresults'] + 'logs/03_multiqc/multiqc.o',
  err = config['pathresults'] + 'logs/03_multiqc/multiqc.e',
  rc = config['pathresults'] + 'logs/03_multiqc/multiqc.rc'
 shell:
  '''
  module load multiqc/1.9
  ( 
   multiqc {input} -o {output}
  ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

############################################## bsmapZ  ####################################
# 		      Step 4 : Alignement bisulfite using bsmapZ			  #
###########################################################################################

rule bsmapz:
 input:
  inp1 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R1_val_1.fq.gz',
  inp2 = config['pathresults'] + '01_trimgalore/{sample}/{sample}_R2_val_2.fq.gz'
 output:
  out = os.path.join(config['pathresults'], '04_bsmapz/{sample}/mapped_{sample}.bam')
 message:
  'mapping reads of reference genome'
 threads: workflow.cores/(config['NBSAMPLE']*1)
 log:
  out = config['pathresults'] + 'logs/04_bsmapz/{sample}/mapping_{sample}.o',
  err = config['pathresults'] + 'logs/04_bsmapz/{sample}/mapping_{sample}.e',
  rc = config['pathresults'] + 'logs/04_bsmapz/{sample}/mapping_{sample}.rc'
 benchmark:
   config['pathresults'] + 'benchmarks/04_bsmapz/{sample}.tsv'
 shell:
  '''
  module load samtools/1.11
  ( 
   {config[tools]}/bsmapz -a {input.inp1} -b {input.inp2} -o {output.out} -d {config[REF_GENOME]} -p {threads}
  ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

###################################### Statistics_Mapping  ####################################
# 				Step 5 : Statistics the mapping				      #
###############################################################################################

rule samtools_stats:
 input:
  inp = os.path.join(config['pathresults'], '04_bsmapz/{sample}/mapped_{sample}.bam')
 output:
  out = os.path.join(config['pathresults'], '05_samtools_stats/{sample}/statistics_{sample}.txt')
 params:
  sta = os.path.join(config['pathresults'] + '05_samtools_stats/{sample}/statistics_{sample}.txt')
 message:
  'Statistics Mapping'
 log:
  out = config['pathresults'] + 'logs/05_samtools_stats/{sample}/stat-alignment_{sample}.o',
  err = config['pathresults'] + 'logs/05_samtools_stats/{sample}/stat-alignment_{sample}.e',
  rc = config['pathresults'] + 'logs/05_samtools_stats/{sample}/stat-alignment_{sample}.rc'
 threads: workflow.cores/(config['NBSAMPLE']*1)
 benchmark:
   config['pathresults'] + 'benchmarks/05_samtools_stats/{sample}.tsv'
 shell:
  '''
  module load samtools/1.11
  (
   samtools stats {input.inp} -r {config[REF_GENOME]} -@ {threads} > {output.out} 
  )1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

################################### Samtools rmdup + Ontarget ###################################
# Step 6 : Elimination of duplicate reads and also to keep only the target regions		#
#################################################################################################

rule samtools_markdup_ontarget:
 input:
  sortbam = os.path.join(config['pathresults'] + '04_bsmapz/{sample}/mapped_{sample}.bam')
 output:
  dupbam = os.path.join(config['pathresults'] + '06_markdup/{sample}/{sample}.dedup.bam'),
  stat = os.path.join(config['pathresults'] + '06_markdup/{sample}/{sample}_statics.txt')
 log:
  out = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.o',
  err = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.e',
  rc = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.rc'
 message:
   'deduplication with samt markdup ontarget'
 threads: workflow.cores/(config['NBSAMPLE']*1)
 benchmark:
  config['pathresults'] + 'benchmarks/06_markdup/{sample}.tsv'
 params:
  t1 = temp(os.path.join(config['pathresults'] + '06_markdup/{sample}/fixmate.bam')),
  t2 = temp(os.path.join(config['pathresults'] + '06_markdup/{sample}/sort.bam')),
  t3 = temp(os.path.join(config['pathresults'] + '06_markdup/{sample}/temp.bam'))
 shell:
  '''
  module load samtools/1.11
  ( 
   samtools fixmate -@ {threads} -O BAM -m {input} {params.t1} && \
   samtools sort -@ {threads} -O BAM {params.t1} -o {params.t2} && \
   samtools markdup -r -@ {threads} -s -f {output.stat} {params.t2} {params.t3} && \
   samtools view -b {params.t3} -L {config[pathbed]} -@ {threads} -h -o {output.dupbam} && \
   samtools index -@ {threads} {output.dupbam} 
  )1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

################################  Samtools rmdup + Not-Ontarget ################################
# Step 6 :  Removing duplicate reads and don't keep target regions 			       #
################################################################################################

rule samtools_markdup:
 input:
  sortbam = os.path.join(config['pathresults'] + '04_bsmapz/{sample}/mapped_{sample}.bam')
 output:
  dupbam = os.path.join(config['pathresults'] + '06_markdup/{sample}/{sample}.dedup.bam'),
  stat = os.path.join(config['pathresults'] + '06_markdup/{sample}/{sample}_statics.txt')
 log:
  out = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.o',
  err = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.e',
  rc = config['pathresults'] + 'logs/06_markdup/{sample}/{sample}_dedup.rc'
 message:
   'deduplication with samt markdup'
 threads: workflow.cores/(config['NBSAMPLE']*1)
 benchmark:
  config['pathresults'] + 'benchmarks/06_markdup/{sample}.tsv'
 params:
  t1 = temp(os.path.join(config['pathresults'] + '06_markdup/{sample}/fixmate.bam')),
  t2 = temp(os.path.join(config['pathresults'] + '06_markdup/{sample}/sort.bam'))
 shell:
  '''
  module load samtools/1.11
  (
   samtools fixmate -@ {threads} -O BAM -m {input} {params.t1} && \
   samtools sort -@ {threads} -O BAM {params.t1} -o {params.t2} && \
   samtools markdup -r -@ {threads} -s -f {output.stat} {params.t2} {output.dupbam} && \
   samtools index -@ {threads} {output.dupbam} 
  )1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

##############################################  methratio.py  ###################################
#   Step 7: Detection of methylated cytosines (mC) in the 3 methylation contexts 		#
#################################################################################################

rule methratio:
 input:
  inp = os.path.join(config['pathresults'] + '06_markdup/{sample}/{sample}.dedup.bam')
 output:
  out = os.path.join(config['pathresults'], '07_methratio/{sample}/meth_{sample}.txt')
 message:
  'methylation caller'
 log:
  out = config['pathresults'] + 'logs/07_methratio/{sample}/{sample}_methratio.o',
  err = config['pathresults'] + 'logs/07_methratio/{sample}/{sample}_methratio.e',
  rc = config['pathresults'] + 'logs/07_methratio/{sample}/{sample}_methratio.rc'
 threads: 
  workflow.cores/(config['NBSAMPLE']*1)
 benchmark:
   config['pathresults'] + 'benchmarks/07_methratio/{sample}.tsv'
 shell:
  '''
   module unload multiqc/1.9
   module unload trimgalore/0.6.5
   module unload fastqc/0.11.9
   module switch python/3.7.6 python/2.7
   module load samtools/1.11
   module load pysam/0.16.0.1
   module list 
  (
   python {config[tools]}/methratio.py {input.inp} -o {output.out} -d {config[REF_GENOME]} -N {threads} -I 
  )1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

################################ split-methylation-context ######################################
# 		step 8 :  Splitting methylation contexts CG,CHG,CHH				#
#################################################################################################

rule splitting_context:
 input:
  inp = os.path.join(config['pathresults'], '07_methratio/{sample}/meth_{sample}.txt')
 output:
  out = os.path.join(config['pathresults'], '08_splitting_context/{sample}/{context}-{sample}.txt')
 log:
  out = config['pathresults'] + 'logs/08_splitting_context/{sample}/{context}-{sample}_splitting.o',
  err = config['pathresults'] + 'logs/08_splitting_context/{sample}/{context}-{sample}_splitting.e',
  rc = config['pathresults'] + 'logs/08_splitting_context/{sample}/{context}-{sample}_splitting.rc'
 benchmark:
   config['pathresults'] + 'benchmarks/08_splitting_context/{context}-{sample}.tsv'
 message:
  'splitting methylation context'
 params:
   var = os.path.join(config['pathresults'], '08_splitting_context/{sample}/{sample}.txt')
 shell:
  '''
  (
   ./script/splitting.sh -i {input.inp} -o {params.var}
  ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
  '''

############################################## preparing_data_R #####################################
# 			Step 9 : Preparing data for statistical analysis			    #
#####################################################################################################

rule preparing_data_R:
 input:
  inp = expand(config['pathresults'] + '08_splitting_context/{sample}/{{context}}-{sample}.txt', sample=SAMPLES)
 output:
  out = expand(config['pathresults'] + '09_preparing_data_R/{{context}}/{{context}}-{sample}.txt', sample=SAMPLES)
 message:
  'preparation context'
 params:
  path = os.path.join(config['pathresults'], '09_preparing_data_R/{context}/')
 shell:
  '''
  mkdir -p {params.path}
  cp {input.inp} {params.path}
  '''

############################################### R methylkit ############################################
# 		Step 10 : Basic statistical analyzes on the detection of methylations		       #
########################################################################################################

rule R_methylkit:
 input:
  expand(config['pathresults'] + '09_preparing_data_R/{{context}}/{{context}}-{sample}.txt', sample=SAMPLES)
 output:
  out = os.path.join(config['pathresults'], '10_methylkit/{context}/{context}.csv'),
  plot = os.path.join(config['pathresults'], '10_methylkit/{context}/{context}.pdf')
 message:
  'stats R in methylkit'
 threads:
  workflow.cores/(config['NBSAMPLE']*3) 
 params:
  path = os.path.join(config['pathresults'], '09_preparing_data_R/{context}/'),
  context = '{context}',
  #coverage = config['COVERAGE'] 
 log: 
  os.path.join(config['pathresults'], 'logs/10_methylkit/{context}.log')
 benchmark:
  config['pathresults'] + 'benchmarks/10_methylkit/{context}.tsv'
 script:
  './script/script_methylkit.R'

########################################################################################################
#					 finally pipeline					       #
########################################################################################################

