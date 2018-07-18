#!/bin/bash

set -e
set -v 
set -x

##Script to call snakefile for bacterial paired-end WGS Illumina data
##Optimized for use on a HPC with SGE scheduler
##aschuerch 062018

##1. Checks
##Check for command line arguments

if [ $# -eq 0 -o "$1" == "-h" -o "$1" == "--help" ]; then
    echo "
###########################################################################
############      Basic microbial WGS analysis pipeline    ################
##                                                                       ##
## for all samples in this folder.                                       ##
##                                                                       ##
## Paired end, compressed sequencing files (fastq.gz)                    ##
## must be present in the same folder from where the script is called.   ##
##                                                                       ##
## Use only the sample names to call the script                          ##
##                                                                       ##
## Example:                                                              ##
##                                                                       ##
## ./bactofidia.sh  ECO-RES-PR1-00001 ECO-RES-PR1-00002                  ##
##                                                                       ##
##                                                                       ##
## or                                                                    ##
##                                                                       ##
## ./bactofidia.sh ALL                                                   ##
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


mkdir -p "$(pwd)"/log
log=$(pwd)/log/call_assembly.txt
touch "$log"
sleep 1

## Check for *fastq.gz

if [ $1 == "ALL" ];then
   files=$(ls *gz | cut -f 1,1 -d _ | uniq | sort -n | tr '\n' ' ')
else
   files="$@"
fi

allfiles=(`echo ${files}`)
echo $allfiles


for i in "$allfiles"
 do
 if [ -f "$i"_*R1*fastq.gz ]
   then
   echo 'Found files for ' "$i" '_*R1*fastq.gz'  2>&1 | tee -a "$log"
  else
   echo 'Sequence files as '"$i"'_*R1*fastq.gz are missing in this folder.
Please execute this script from the location of the sequencing files or exclude 
the sample.
Exiting.' 2>&1 | tee -a "$log"
   exit 0
  fi
 done


# check if conda is installed
if command -v conda > /dev/null; then
 echo  2>&1| tee -a "$log"
else
 echo "Miniconda missing. Installing...." 
 wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
 chmod +x Miniconda3-latest-Linux-x86_64.sh
 mkdir -p ~/tmp
 ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/tmp/Miniconda3
 rm Miniconda3-latest-Linux-x86_64.sh
 export PATH=~/tmp/Miniconda3/bin:$PATH 
 export PYTHONPATH=~/tmp/Miniconda3/pkgs/
 conda config --add channels conda-forge
 conda config --add channels defaults
 conda config --add channels r
 conda config --add channels bioconda
 export PERL5LIB="~/tmp/Miniconda3/lib/perl5/site_perl/5.22.0"
fi

# Check and activate snakemake 
# source activate snakemake || echo "Please create a virtual environment with snakemake and python3 with 'conda create -n snakemake snakemake python=3.5"

echo |  2>&1 tee -a "$log"
echo "The results will be generated in this location: " 2>&1| tee -a "$log"
echo "$(pwd)"/results 2>&1| tee -a "$log"
echo |  2>&1 tee -a "$log"
sleep 1

echo "The logfiles will be generated here: " 2>&1 | tee -a "$log"
echo "$(pwd)"/log  2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"
sleep 1

# determine read length and config files

for i in "${files%% *}"
  do
    length=$(zcat "$i"_*R1*fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort | uniq -c | sort -rn | head -n 1 | rev | cut -f 1,1 -d " "| rev)
  done

if [[ "$length" == 151 ]];then
  configfile=config.yaml
elif [[ "$length" == 251 ]]; then
  configfile=config_miseq.yaml
else
  echo 'please provide a custom config file (e.g. config_custom.yaml): '
  read -r configfile
fi


# Check virtual environment or create
#virtenv=$(grep -A 1 virtual_environment config.yaml | grep -v "#" | cut -f 2,2 -d ":" | sed -e 's/^[[:space:]]*//')
#echo "$virtenv"
#conda env list | grep "$virtenv" || conda create -y -n "$virtenv" --file package-list.txt


# Check and activate snakemake or create
source activate snakemake || conda create -y -n snakemake snakemake


echo 2>&1 |tee -a "$log"
echo "Read length was determined as: " 2>&1| tee -a "$log"
echo "$length" 2>&1| tee -a "$log"
echo "$configfile" "will be used as configfile"   2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"


sleep 1
 
# concatenate for rev and put into data/ folder:
mkdir -p data

while IFS=' ' read -ra samples
 do
   for i in "${samples[@]}";
    do
     echo "$i"
     cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
     cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz
  done
 done <<< "$files"


#check if it is on hpc

if command -v qstat > /dev/null; then

#get e-mail to send the confirmation to
 emaildict=/hpc/dla_mm/data/shared_data/bactofidia_config/email.txt
 if [[ -e "$emaildict" ]]; then
   echo 'Email file found' 2>&1| tee -a "$log"
   while read name mail
    do
      if [[ "$name" == "$(whoami)" ]]; then
       email="$mail"
      fi 
    done < "$emaildict"
 else
   echo 'please provide your e-mail '
   read -p email
 fi

echo 'An e-mail will be sent to '"$email"' upon job completion.' 2>&1| tee -a "$log" 

#command on cluster (SGE)
 snakemake \
 --snakefile Snakefile.assembly \
 --latency-wait 60 \
 --config configfile="$configfile" \
 --verbose \
 --printshellcmds \ 
 --forceall \
 --keep-going \
 --restart-times 5\
 --use-conda \
 --cluster \
 'qsub -V -cwd -l h_vmem=125G -l h_rt=04:00:00 -e log/ -o log/ ' \
 --jobs 100 2>&1| tee -a "$log"

#job to send an e-mail
job=log/bactofidia_done.sh
touch "$job"
echo "#!/bin/bash" > "$job"
echo "sleep 1" > "$job"

echo qsub -e "$(pwd)"/log/ -o "$(pwd)"/log/ -m ae -M "$email" "$job"
qsub -e "$(pwd)"/log/ -o "$(pwd)"/log/ -m ae -M "$email" "$job"

else

#if not on a cluster
echo "snakemake --snakefile Snakefile.assembly --keep-going --config configfile=""$configfile"
snakemake --snakefile Snakefile.assembly --keep-going --config configfile="$configfile"

#for the CI
if [ $? -eq 0 ]
then
  echo "Successfully finished job"
  exit 0
else
  echo "Could not finish job" >&2
  exit 1
fi

fi

