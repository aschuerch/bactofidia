#!/bin/bash

# The path to all the binaries
PATH1=/hpc/local/CentOS7/dla_mm/bin/

##re-written and simplified version, Anita june 2016

##Checks
##Check for command line arguments

if [ $# -eq 0 ]; then
    echo "
###########################################################################
## #######      Assembly pipeline       ###########                      ##
## for samples in the same folder. The sample names should be given      ##
## as command line options. Compressed sequencing files (.gz)            ##
## must be present in the same folder from where the script is called.   ##
##                                                                       ##
## Example                                                               ##
##                                                                       ##
## call_assembly.sh ECO-RES-PR1-00001 ECO-RES-PR1-00002 ...              ##
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
## Anita Schurch May 2016                                                ##
###########################################################################"
    exit 
fi


## path to the results
SHARED_OUTPUT=/hpc/dla_mm/automatedassembly
CURRENTDIR=$(pwd)
timestamp=$(date +%s)

## copy AssemblyParameters if requested.
if [[ $1 == "AssemblyParameters" ]]; then
 cp $SHARED_OUTPUT/parameterfiles/AssemblyParameter.txt $CURRENTDIR/AssemblyParameter_$timestamp.txt
 exit 0
fi 




SCRIPTDIR=/hpc/shared/dla_mm/Assembly_pipeline
#### SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
user=$(whoami)

declare -A emaildict=( ["aschurch"]="a.c.schurch@umcutrecht.nl" \
["jbayjanov"]="J.Bayjanov@umcutrecht.nl" \
["mrogers"]="M.R.C.Rogers-2@umcutrecht.nl" \
["salonso"]="S.ArredondoAlonso@umcutrecht.nl")
email="${emaildict[$user]}"


##Check for *fastq.gz
for i in "$@"
 do
  ls -1 "$CURRENTDIR"/"$i"_*R1*fastq.gz > /dev/null 2>&1
  if [[ "$?" = "0" ]]; then
   echo "Found file(s) for "$i
  else
   echo "Sequence file(s) as    $i*fastq.gz          is/are missing in this folder.
Please execute this script from the location of the sequencing files or exclude the sample.
Exiting."
   exit 0
  fi
 done

##use first sample as project name

dir=P-$1
dir=$dir"_"$timestamp


log=$SHARED_OUTPUT/logs/"$dir"_"$user"_call_assembly.txt
mkdir -p $SHARED_OUTPUT/$dir
touch $log
sleep 1

echo |  2>&1 tee -a $log
echo "The results will be generated in this location: " 2>&1| tee -a $log
echo $SHARED_OUTPUT/$dir  2>&1| tee -a $log
echo |  2>&1 tee -a $log
sleep 1

echo "The logfiles will be generated here: " 2>&1 | tee -a $log
echo $log  2>&1| tee -a $log
echo $SHARED_OUTPUT/logs/"$dir"_"$user"_samplepipe.txt   2>&1| tee -a $log
echo $SHARED_OUTPUT/logs/"$dir"_"$user"_projectpipe.txt    2>&1| tee -a $log
echo 2>&1 |tee -a $log
sleep 1

if [[ -z $email ]];then
echo  2>&1| tee -a $log
echo "User not found. Please enter the e-mail address to which the confirmation of completion needs to be send"
read email
fi
echo  2>&1| tee -a $log
echo "Confirmation of completion will be sent to: "$email  2>&1| tee -a $log
echo  2>&1| tee -a $log
sleep 1


genusfile=$SCRIPTDIR/parameterfiles/genusfile.txt
findgenus () {
	grep "^$1:" $genusfile | sed s/"$1:"//
}
genus=$(findgenus $(echo "$1" |cut -f 1,1 -d "-"))

if [[ -z $genus ]];then
  echo "Genus will be determined from parameterfile"   2>&1| tee -a $log
  echo  2>&1| tee -a $log
fi


echo "1- Copying sequencing files"  2>&1| tee -a $log
echo  2>&1| tee -a $log
sleep 1
##symlink data to $SHARED_OUTPUT/$dir
for i in "$@" 
do
ln -s $CURRENTDIR/$i*fastq.gz $SHARED_OUTPUT/$dir/
done

echo "2- Creating parameter file"  2>&1| tee -a $log
echo  2>&1| tee -a $log
#create parameterfile
paramfile=$SHARED_OUTPUT/$dir/AssemblyParameter_"$timestamp".txt
if [[ -e  $paramfile ]];
then
   rm $paramfile
fi

touch $paramfile

sleep 1

##Populate parameter file
echo " " >> $paramfile
echo "genus:"$genus >> $paramfile
echo "user:"$user >> $paramfile
echo "email:"$email >> $paramfile
echo "timestamp:"$timestamp  >> $paramfile
echo "project:"$dir >> $paramfile
echo "samples:"$@ >> $paramfile
echo " " >> $paramfile
 
ls -1 "$CURRENTDIR"/AssemblyParameter*txt > /dev/null 2>&1
  if [[ "$?" = "0" ]]; 
then
    echo "AssemblyParameter file found. Custom parameters will be used." 2>&1| tee -a $log
    cat $CURRENTDIR/AssemblyParameter*.txt  >> $paramfile
else 
    echo "No or several new AssemblyParameter.txt files found. Standard parameters used for pipeline."  2>&1| tee -a $log
    cat $SCRIPTDIR/parameterfiles/AssemblyParameter.txt >> $paramfile
fi

sleep 1
echo  2>&1| tee -a $log
echo "The parameterfile can be found here: "$paramfile 2>&1|tee -a $log

echo  2>&1| tee -a $log
echo "3- Started job creation" 2>&1 |tee -a $log
echo  2>&1| tee -a $log
echo "The number of jobs for assembly is: "$#  2>&1| tee -a $log
echo  2>&1| tee -a $log


echo "qsub -t 1-$# -N J$timestamp -o $SHARED_OUTPUT/logs/"$dir"_"$user"_samplepipe.txt -e $SHARED_OUTPUT/logs/"$dir"_"$user"_samplepipe.txt $SCRIPTDIR/samplepipe.sh "$paramfile" " >> $log

##array job for assembly pipeline
qsub -t 1-$# -N J$timestamp \
-o $SHARED_OUTPUT/logs/"$dir"_"$user"_samplepipe.txt \
-e $SHARED_OUTPUT/logs/"$dir"_"$user"_samplepipe.txt \
$SCRIPTDIR/samplepipeTEST.sh "$paramfile"  
 

echo "Jobnumber of the array job is: " J$timestamp 2>&1| tee -a $log


echo "qsub -N Fin_"$dir" -e  $SHARED_OUTPUT/logs/roundup_"$dir"_"$user" -o $SHARED_OUTPUT/logs/roundup_"$dir"_"$user" -M $email -m e -hold_jid J$timestamp $SCRIPTDIR/projectpipe.sh "$paramfile" " >> $log

##roundup job
qsub -N Fin_"$dir" \
-e  $SHARED_OUTPUT/logs/"$dir"_"$user"_projectpipe.txt \
-o $SHARED_OUTPUT/logs/"$dir"_"$user"_projectpipe.txt \
-M $email \
-m e \
-hold_jid J$timestamp \
$SCRIPTDIR/projectpipe.sh "$paramfile"
