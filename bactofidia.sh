#!/bin/bash

set -e

##Script to call snakefile for bacterial paired-end WGS Illumina data
##aschuerch 052017

##1. Checks
##Check for command line arguments

if [ $# -eq 0 ]; then
    echo "
###########################################################################
############      Basic microbial WGS analysis pipeline    ################
## for all available samples in this folder.                             ##
## Compressed sequencing files (.gz)                                     ##
## must be present in the same folder from where the script is called.   ##
## Use only the sample name to call the script                           ##
##                                                                       ##
## Example                                                               ##
##                                                                       ##
## bactofidia.sh  ECO-RES-PR1-00001 ECO-RES-PR1-00002                    ##
##                                                                       ##
##                                                                       ##
## Before running the pipeline for the first time, a virtual             ##
## environment needs to be created. Packages and versions are specified  ##
## in package-list.txt. Adjust this file to your needs. See              ##
## bioconda.github.io for available packages.                            ##
##                                                                       ##
## Create the environment with                                           ##
##                                                                       ##
## conda create --file package-list.txt -n [bactofidia_custom]           ##
##                                                                       ##
## where [bactiofidia_custom] matches the name of the environment        ##
## given in the config.yaml. The config.yaml file can be adjusted for    ##
## parameters of the different tools.If parameters different from the    ##
## standard parameters will be used,                                     ##
## you can adjust package-list.txt config.yaml                           ##
## to your needs before running this script                              ##
##                                                                       ##
##                                                                       ##
## Anita Schurch May 2017                                                ##
###########################################################################"
    exit
fi


## Check for *fastq.gz
for i in "$@"
 do
  if (ls -1 "$i"_*R1*fastq.gz > /dev/null 2>&1)
   then
   echo 'Found file(s) for ' "$i"
  else
   echo 'Sequence file(s) as '"$i"'*fastq.gz          is/are missing in this folder.
Please execute this script from the location of the sequencing files or exclude 
the sample.
Exiting.'
   exit 0
  fi
 done

mkdir -p $(pwd)/log
log=$(pwd)/log/call_assembly.txt
touch "$log"
sleep 1


# check if conda is installed
if command -v conda > /dev/null; then
 echo  2>&1| tee -a "$log"
else
 echo "Miniconda missing" 
 exit 0
fi


# Check and activate snakemake 
source activate snakemake || echo "Please create a virtual environment with snakemake and python3" 


echo |  2>&1 tee -a "$log"
echo "The results will be generated in this location: " 2>&1| tee -a "$log"
pwd 2>&1| tee -a "$log"
echo |  2>&1 tee -a "$log"
sleep 1

echo "The logfiles will be generated here: " 2>&1 | tee -a "$log"
echo "$(pwd)"/log  2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"
sleep 1


# determine if miseq or hiseq
for i in "$@"
 do
 num=$(find . -maxdepth 1 -name "$i" | wc -l)
 if [[ $((num/2)) -gt 1 ]];then
  echo  2>&1| tee -a "$log"
  echo "Sequencing type is NextSeq"  2>&1| tee -a "$log"
  seq='nextseq'
else
  echo  2>&1| tee -a "$log"
  echo "Sequencing type is MiSeq"  2>&1| tee -a "$log"
  seq='miseq'
fi
done


# concatenate for rev and put into data/ folder:
mkdir -p data
for i in "$@"
 do
   cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
   cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz
done


#check if it is on hpc
if command -v qstat > /dev/null; then

 if [[ $seq == 'nextseq' ]]; then

 echo "snakemake \
 --latency_wait 60 \
 --config krange="33,55,71" \
 --verbose \
 --forceall \
 --cluster \
 --keepgoing \
 --restart_times 5\
 'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
 --jobs 100 "

 else #miseq

 echo "snakemake \
 --latency_wait 60 \
 --config krange="57,97,127" \
 --verbose \
 --keepgoing \
 --restart_times 5\
 --cluster \
 --forceall \
 'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
 --jobs 100 " 

 snakemake \
 --latency_wait 60 \
 --config krange="57,97,127" \
 --verbose \
 --forceall \
 --keepgoing \
 --restart_times 5\
 --cluster \
 'qsub -cwd -l h_vmem=125G -l h_rt=04:00:00 -e log/ -o log/ -M a.c.schurch@umcutrecht.nl ' \
 --jobs 100

 fi
else

snakemake --keep-going --forceall


fi

