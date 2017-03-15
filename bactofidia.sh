#!/bin/bash

##Script to call snakefile for bacterial paired-end WGS Illumina data

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
## If parameters different from the standard parameters will be used,    ##
## you can adjust config.yaml to your needs before running this script   ##
##                                                                       ##
 #                                                                       ##
## Anita Schurch Feb 2017                                                ##
###########################################################################"
    exit
fi


##Check for *fastq.gz
for i in "$@"
 do
  ls -1 "$i"_*R1*fastq.gz > /dev/null 2>&1
  if [[ "$?" = "0" ]]; then
   echo "Found file(s) for "$i
  else
   echo "Sequence file(s) as    $i*fastq.gz          is/are missing in this folder.
Please execute this script from the location of the sequencing files or exclude 
the sample.
Exiting."
   exit 0
  fi
 done

#if needed download miniconda installation
if command -v conda > /dev/null; then
 echo  2>&1| tee -a $log

else

 wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
 chmod +x Miniconda2-latest-Linux-x86_64.sh
 mkdir -p ~/tmp
 ./Miniconda2-latest-Linux-x86_64.sh -b -p ~/tmp/Miniconda2
 rm Miniconda2-latest-Linux-x86_64.sh
 export PATH=~/tmp/Miniconda2/bin:$PATH 
 export PYTHONPATH=~/tmp/Miniconda2/pkgs/
 conda config --add channels conda-forge
 conda config --add channels defaults
 conda config --add channels r
 conda config --add channels bioconda
 conda create -y -n snakemake python=3.5 snakemake

fi


#activate snakemake env
source activate snakemake


log=$(pwd)/log/call_assembly.txt
mkdir -p log
touch $log
sleep 1

echo |  2>&1 tee -a $log
echo "The results will be generated in this location: " 2>&1| tee -a $log
echo $(pwd)  2>&1| tee -a $log
echo |  2>&1 tee -a $log
sleep 1

echo "The logfiles will be generated here: " 2>&1 | tee -a $log
echo $(pwd)/log  2>&1| tee -a $log
echo 2>&1 |tee -a $log
sleep 1

#if [[ -z $email ]];then
#echo  2>&1| tee -a $log
#echo "User not found. Please enter the e-mail address to which the confirmation of completion needs to be send
#"
#read email
#fi


# determine if miseq or hiseq
seq='miseq'
if [[ 2*$# > $(ls -1 *fastq.gz) ]];then
  echo  2>&1| tee -a $log
  echo "Sequencing type is NextSeq"  2>&1| tee -a $log

  seq='nextseq'
else
  echo  2>&1| tee -a $log
  echo "Sequencing type is MiSeq"  2>&1| tee -a $log

  seq='miseq'
fi



# concatenate for rev and put into data/ folder:
mkdir -p data
for i in "$@"
 do
   cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
   cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz
done


#check if it is on hpc
if command -v qstat > /dev/null; then

 if [[ seq == 'nextseq' ]]; then

 echo "snakemake \
 --latency-wait 60 \
 --config krange="33,55,71" \
 --verbose \
 --cluster \
 'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
 --jobs 100 "

 else #miseq

 echo "snakemake \
 --latency-wait 60 \
 --config krange="57,97,127" \
 --verbose \
 --keep-going \
 --cluster \
 'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
 --jobs 100 " 

 snakemake \
 --latency-wait 60 \
 --config krange="57,97,127" \
 --verbose \
 --keep-going \
 --cluster \
 'qsub -cwd -l h_vmem=125G -l h_rt=04:00:00 -e log/ -o log/ -M a.c.schurch@umcutrecht.nl ' \
 --jobs 100

 fi
else

snakemake --keep-going 

fi

