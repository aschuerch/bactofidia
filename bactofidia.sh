#!/bin/bash -i

##to debug
#set -e
#set -v
#set -x

########################################
##Script to call snakefile for bacterial paired-end WGS Illumina data
##Run on a HPC with SLURM scheduler by requesting a node with sufficient memory 
## screen -S [session_name]
## srun --time=24:00:00 --mem=32G --pty bash --output=log/ --error log/
##aschuerch 042020
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
## Version March2020                                                     ##
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

#Check if conda is installed, if not found attempt to install in a temporary folder

if command -v conda > /dev/null; then
 echo  2>&1| tee -a "$log"
 echo 
 echo "conda found" | tee -a "$log"
else
 echo
 echo "conda missing"
 echo "Install Miniconda with" 
 echo
 echo "    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
 echo "    chmod +x Miniconda3-latest-Linux-x86_64.sh"
 echo "    ./Miniconda3-latest-Linux-x86_64.sh"
 echo "and follow the prompts."
 echo "After installation, configure the channels with"
 echo
 echo "    conda config --add channels defaults"
 echo "    conda config --add channels bioconda"
 echo "    conda config --add channels conda-forge"
 exit 1
fi

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

# Write to log
echo 2>&1 |tee -a "$log"
echo "Read length was determined as: " 2>&1| tee -a "$log"
echo "$length" 2>&1| tee -a "$log"
echo "$configfile" "will be used as configfile"   2>&1| tee -a "$log"
echo 2>&1 |tee -a "$log"

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo 2>&1 |tee -a "$log"
echo "snakemake found" 2>&1 |tee -a "$log"
else
echo 2>&1 |tee -a "$log"
echo "snakemake will be installed" 2>&1 |tee -a "$log"
conda install -y snakemake
fi


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

#if command -v qstat > /dev/null; then

#get e-mail to send the confirmation to
# emaildict=/hpc/dla_mm/data/shared_data/bactofidia_config/email.txt
# if [[ -e "$emaildict" ]]; then
#   echo 'Email file found' 2>&1| tee -a "$log"
#   while read name mail
 #   do
  #    if [[ "$name" == "$(whoami)" ]]; then
  #     email="$mail"
  #    fi 
  #  done < "$emaildict"
# else
#   echo 'please provide your e-mail '
#   read -p email
# fi

#echo 'An e-mail will be sent to '"$email"' upon job completion.' 2>&1| tee -a "$log" 

#command on cluster (SLURM)
# snakemake \
# --snakefile Snakefile.assembly \
# --profile config/slurm
# --config configfile="$configfile" \
# --verbose \
# --forceall \
# --keep-going \
# --restart-times 3\
# --use-conda \
# --jobs 10 2>&1| tee -a "$log"

#job to send an e-mail
#job=log/bactofidia_done.sh
#touch "$job"
#echo "#!/bin/bash" > "$job"
#echo "sleep 1" > "$job"

#echo qsub -e "$(pwd)"/log/ -o "$(pwd)"/log/ -m ae -M "$email" "$job"
#qsub -e "$(pwd)"/log/ -o "$(pwd)"/log/ -m ae -M "$email" "$job"

#else

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

#fi
