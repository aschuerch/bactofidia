#!/bin/bash

#calls snakefile for bacterial paired-end WGS Illumina data


#activate snakemake env
source activate snakemake-tutorial

##Checks
##Check for command line arguments

if [ $# -eq 0 ]; then
    echo "
###########################################################################
## #######      Assembly pipeline       ###########                      ##
## for all available samples in this folder.                             ##
## Compressed sequencing files (.gz)                                     ##
## must be present in the same folder from where the script is called.   ##
##                                                                       ##
## Example                                                               ##
##                                                                       ##
## call_assembly.sh  ECO-RES-PR1-00001 ECO-RES-PR1-00002                 ##
##                                                                       ##
##                                                                       ##
## If parameters different from the standard parameters will be used,    ##
## you can copy AssemblyParameters.txt file by calling the script with   ##
##                                                                       ##  
## call_assembly.sh AssemblyParameters                                   ##
##                                                                       ##
## adjust it to your needs,save it as AssemblyParameters_[version].txt   ##
## Place it in the same folder in which your sequencing files are.       ##
##                                                                       ##     
              
## Anita Schurch Feb 2017                                                ##
###########################################################################"
    exit 
fi


## copy AssemblyParameters if requested.
#if [[ $1 == "AssemblyParameters" ]]; then
# cp $SHARED_OUTPUT/parameterfiles/AssemblyParameter.txt $CURRENTDIR/AssemblyParameter_$timestamp.txt
# exit 0
#fi 

#defines user, determines who to send the email
#user=$(whoami)

#declare -A emaildict=( ["aschurch"]="a.c.schurch@umcutrecht.nl" \
#["jbayjanov"]="J.Bayjanov@umcutrecht.nl" \
#["mrogers"]="M.R.C.Rogers-2@umcutrecht.nl" \
#["salonso"]="S.ArredondoAlonso@umcutrecht.nl")
#email="${emaildict[$user]}"


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
  echo "Sequencing type is HiSeq"
  seq='hiseq'
else
  echo  2>&1| tee -a $log
  echo "Sequencing type is MiSeq"
  seq='miseq'
fi



# concatenate for rev and put into data/ folder:
mkdir -p data
for i in "$@"
 do
   cat "$i"*R1*.fastq.gz > data/"$i"_R1.fastq.gz
   cat "$i"*R2*.fastq.gz > data/"$i"_R2.fastq.gz
done

# create yaml 
# start snakefile

if [[ seq == 'hiseq' ]]; then

echo "snakemake \
--latency-wait 60 \
--config krange="33,55,71" \
--verbose \
--cluster \
'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
--jobs 100 "

snakemake -npr

else

echo "snakemake \
--latency-wait 60 \
--config krange="57,97,127" \
--verbose \
--cluster \
'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
--jobs 100 " 

snakemake \
--latency-wait 60 \
--config krange="57,97,127" \
--verbose \
--cluster \
'qsub -cwd -l h_vmem=48G -l h_rt=04:00:00 -e log/ -o log/ ' \
--jobs 100

snakemake -npr
fi
