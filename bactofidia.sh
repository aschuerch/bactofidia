#!/bin/bash -i

##to debug
#set -e
#set -v 
#set -x

########################################
##Script to call snakefile for bacterial paired-end WGS Illumina data
##Optimized for use on a HPC with SGE scheduler
##aschuerch 082018
########################################


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
##                                                                       ##
## Example:                                                              ##
##                                                                       ##
## ./bactofidia.sh  ECO-RES-PR1-00001_R1.fastq.gz ECO-RES-PR1-00001_R2.fastq.gz
##  ECO-RES-PR1-00002_R1.fastq.gz ECO-RES-PR1-00002_R2.fastq.gz          ##
##                                                                       ##
## or                                                                    ##
##                                                                       ##
## ./bactofidia.sh ALL                                                   ##
##                                                                       ##
## Packages and versions are specified in envs/packages.yml.             ## 
## See bioconda.github.io for available packages.                        ##
## Command line parameters for individual tools can be adjusted in       ##
## config/config.yaml                                                    ##
##                                                                       ##
## Version Aug2018                                                       ##
###########################################################################"
    exit
fi

## create logs

mkdir -p "$(pwd)"/log
log=$(pwd)/log/call_assembly.txt
touch "$log"
sleep 1

## Check for *fastq.gz files



if [ "$1" == """ALL""" ];then
   echo "All fastq.gz files will be processed"  2>&1 | tee -a "$log"
   files=(./*fastq.gz)
else
   echo else
   files=( "$@" )
fi


for file in "${files[@]}"
do
  if [ -e "$file" ]
   then   # Check whether file exists.
     echo 'Found files for ' "$file"  2>&1 | tee -a "$log"
   else
     echo 'Sequence files as '"$file"'_*R1*fastq.gz are missing in this folder.
  Please execute this script from the location of the sequencing files or exclude the sample.
 Exiting.' 2>&1 | tee -a "$log"
 exit 1
   fi
done

# check if results directory already exist, wait for user input to continue
if [ -d results ] 
  then 
    read -p "The results directory already exists and will be overwritten. Do you want to continue (Y/N)? " -n 1 -r 
    echo    # new line
    if [[ $REPLY =~ ^[Nn]$ ]]
     then
      exit 1
    else
     rm -r results/ 2>&1 | tee -a "$log"
     if [ -d data ]
     then
       rm -r data/ 2>&1 | tee -a "$log"
     fi
    fi
fi


# check if conda is installed, if not found attempt to install in a temporary folder

if command -v conda > /dev/null; then
 echo  2>&1| tee -a "$log"
else
 echo "Miniconda missing. Installing...." 
 wget https://repo.continuum.io/miniconda/Miniconda3-4.3.31-Linux-x86_64.sh
 chmod +x Miniconda3-latest-Linux-x86_64.sh
 mkdir -p ~/tmp
 ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/tmp/Miniconda3
 rm Miniconda3-latest-Linux-x86_64.sh
 export PATH=~/tmp/Miniconda3/bin:$PATH 
 export PYTHONPATH=~/tmp/Miniconda3/pkgs/
 conda config --add channels defaults
 conda config --add channels bioconda
 conda config --add channels conda-forge
 export PERL5LIB=~/tmp/Miniconda3/lib/perl5/site_perl/5.22.0
fi

echo |  2>&1 tee -a "$log"
echo "The results will be generated in this location: " 2>&1| tee -a "$log"
echo "$(pwd)"/results 2>&1| tee -a "$log"
echo |  2>&1 tee -a "$log"
sleep 1

echo "The logfiles will be generated here: " 2>&1 | tee -a "$log"
echo "$(pwd)"/log  2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"
sleep 1

# determine read length

length=$(zcat "${files[0]}" | awk '{if(NR%4==2) print length($1)}' | sort | uniq -c | sort -rn | head -n 1 | rev | cut -f 1,1 -d " "| rev)

# determine which config file according to read length (only 251 and 151 implemented a the moment)

if [[ "$length" == 151 ]];then
  configfile=config/config.yaml
elif [[ "$length" == 251 ]]; then
  configfile=config/config_miseq.yaml
else
  echo 'Sequence length is '"$length"', please provide a custom config file (e.g. config_custom.yaml): '
  read -r configfile
fi


# Check and activate snakemake or create new environment
source activate snakemake || conda activate snakemake || conda create -y -n snakemake snakemake=5.2.2 python=3.5 && source activate snakemake || conda activate snakemake
# Write to log

echo 2>&1 |tee -a "$log"
echo "Read length was determined as: " 2>&1| tee -a "$log"
echo "$length" 2>&1| tee -a "$log"
echo "$configfile" "will be used as configfile"   2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"


sleep 1
 
# concatenate forward and reverse and put into data/ folder:

mkdir -p data

for file in "${files[@]}"
 do
  i="${file%%_*}"
  echo "$i"
  cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
  cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz 
 done

#check if it is on hpc (with sge scheduler)

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
 --forceall \
 --keep-going \
 --restart-times 5\
 --use-conda \
 --cluster \
 'qsub -V -cwd -l h_vmem=32G -l h_rt=04:00:00 -e log/ -o log/ ' \
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
echo "snakemake --snakefile Snakefile.assembly --use-conda --printshellcmds --keep-going --config configfile=""$configfile"
snakemake --snakefile Snakefile.assembly --use-conda --printshellcmds  --keep-going --config configfile="$configfile"

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
