#!/bin/bash
#activate snakemake env
source activate snakemake-tutorial
# concatenate for rev and put into data/ folder:
# cat *all lanesR1.fastq.gz > R1.fastq.gz
# create yaml 
# start snakefile

#snakemake -np
snakemake --verbose --cluster 'qsub -cwd -l h_vmem=24G -l h_rt=04:00:00 -e ~/logs/ -o ~/logs/' --jobs 15
