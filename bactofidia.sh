#!/bin/bash
set -e

##Script to call snakefile for bacterial paired-end WGS Illumina data
##aschuerch 092017

##1. Checks
##Check for command line arguments

if [ $# -eq 0 ]; then
    echo "
###########################################################################
############      Basic microbial WGS analysis pipeline    ################
##                                                                       ##
## for all samples in this folder.                                       ##
##                                                                       ##
## Compressed sequencing files (fastq.gz)                                ##
## must be present in the same folder from where the script is called.   ##
##                                                                       ##
## Use only the sample names to call the script                          ##
##                                                                       ##
## Example:                                                              ##
##                                                                       ##
## bactofidia.sh  ECO-RES-PR1-00001 ECO-RES-PR1-00002                    ##
##                                                                       ##
##                                                                       ##
## Before running the pipeline for the first time, a virtual             ##
## environment needs to be created. Packages and versions are specified  ##
## in package-list.txt. See bioconda.github.io for available packages.   ##
##                                                                       ##
## Create the environment with                                           ##
##                                                                       ##
## conda create --file package-list.txt -n bactofidia_standard201709     ##
##                                                                       ##
##                                                                       ##
## Anita Schurch Aug 2017                                                ##
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

mkdir -p "$(pwd)"/log
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

# determine read length and config file

for i in "$@"
 do
 length=$(zcat "$i"_*R1*fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort | uniq -c | sort -rn | head -n 1 | rev | cut -f 1,1 -d " "| rev)
 done

if [[ "$length" == 151 ]];then
  configfile=config.yaml
elif [[ "$length" == 251 ]]; then
  configfile=config_miseq.yaml
else
  echo 'please provide a custom config file (e.g. config_custom.yaml) '
  read -r configfile
fi

echo 2>&1 |tee -a "$log"
echo "Read length was determined as: " 2>&1| tee -a "$log"
echo "$length" 2>&1| tee -a "$log"
echo "$configfile" "will be used as configfile"   2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"


sleep 1
 
# concatenate for rev and put into data/ folder:
mkdir -p data
for i in "$@"
 do
   cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
   cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz
done


#check if it is on hpc
if command -v qstat > /dev/null; then

 snakemake \
 --latency-wait 60 \
 --config configfile="$configfile" \
 --verbose \
 --forceall \
 --keep-going \
 --restart-times 5\
 --cluster \
 'qsub -cwd -l h_vmem=125G -l h_rt=04:00:00 -e log/ -o log/ -M a.c.schurch@umcutrecht.nl ' \
 --jobs 100 2>&1| tee -a "$log"

else

snakemake --keep-going --config configfile="$configfile" 2>&1| tee -a "$log"


fi

