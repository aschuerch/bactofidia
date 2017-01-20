#!/bin/bash

#set -x
#set -u

#$ -cwd
#$ -m sa             # When to send alerts b=beginning, e=end, a=abort, s=suspend
#$ -M a.c.schurch@umcutrecht.nl # email address for alerts
#$ -l h_vmem=24G
#$ -l h_rt=06:00:00

SHARED_OUTPUT=/hpc/dla_mm/automatedassembly

echo task_id:$SGE_TASK_ID

sleep $((60*$SGE_TASK_ID))

#activate virtual env
export CONDA_OLD_PS1="[\u@\h \W]\$";

paramfile=$1
parameter () {
	grep "^$1:" $paramfile | sed s/"$1:"//
}

samples=$(parameter samples)
IFS=' ' read -ra ADDR <<< $samples
index=$(expr $SGE_TASK_ID - 1)
samplename="${ADDR[$index]}"

echo samplename:$samplename
user=$(parameter user)
echo user:$user
project=$(parameter project)
echo project:$project
##repository were all the files are stored
DIR=$SHARED_OUTPUT/$project


genus=$(parameter genus)
if [[ -z $genus ]];then
  echo "Genus is empty"
else
  echo "Genus:"$genus
fi

#version of the assembly pipeline
V=$project
echo Version:$V

## folder where analyses takes place
assemblyfolder=$SHARED_OUTPUT/assemblies/$project/"$samplename"
mkdir -p $SHARED_OUTPUT/assemblies/$project
echo "assemblyfolder:"$assemblyfolder
mkdir -p $assemblyfolder

## directory to put all results
resultsdir=$DIR/$project
mkdir -p $resultsdir
echo "resultsfolder:"$resultsdir
cp $paramfile $resultsdir/

#minimum length of contigs
echo "minimum length: $(parameter minlen)
#minimum coverage depth of contigs
echo "minimum coverage: $(parameter mincov)
echo ""


#########Script starts here#############

echo ""
for r in R1 R2
do
echo "concatenate forward and reverse lanes for " $samplename.$r
if [[ ! -f $DIR/$samplename.$r.untrimmed.fastq.gz ]]; then
 cat $DIR/"$samplename"*$r*gz > $DIR/$samplename.$r.untrimmed.fastq.gz
fi
done


echo ""
echo trimming:$(parameter trimming)
if [[ $(parameter trimming) = no ]]; then
for r in R1 R2
do
echo "no trimming"
seqtk seq $DIR/$samplename.$r.untrimmed.fastq.gz > $DIR/$samplename.$r.fastq
done

else
echo "trimming"
trimstat="$resultsdir"/Trimmingstats_$project.tsv

if [[ ! -f $trimstat ]]; then
touch $trimstat
echo -e "Reads\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high" >> $trimstat
fi
for r in R1 R2
do
trimstat="$resultsdir"/Trimmingstats_$project.tsv
echo "trimming with seqtk trimfq $(parameter seqtkparam) $DIR/$samplename.$r.untrimmed.fastq.gz > $DIR/$samplename.$r.fastq"
seqtk trimfq "$(parameter seqtkparam)" $DIR/$samplename.$r.untrimmed.fastq.gz > $DIR/$samplename.$r.fastq
echo -n ""$samplename"-"$r" untrimmed" >> $trimstat
seqtk fqchk $DIR/$samplename.$r.untrimmed.fastq.gz | grep ALL  | sed s/ALL//g >> $trimstat

echo -n ""$samplename"-"$r" trimmed"  >> $trimstat
seqtk fqchk $DIR/$samplename.$r.fastq | grep ALL   | sed s/ALL//g >> $trimstat
done
fi


#echo "remove untrimmed files to save space"
#rm $DIR/*gz


echo ""
if [[ $(parameter assembly) = yes ]]; then
echo "assembly: "
## determine read length with seqtk
LEN=$(seqtk fqchk "$DIR"/"$samplename".R1.fastq | grep max_len | cut -f 2 -d ";" | cut -f 2 -d ":")
echo "read length:" $LEN
## determine kmer range with python script
krange=''
if [[ $(parameter krange) = default ]]; then
   if [[  $LEN > 151 ]]; then
       krange=$krange'57,97,127'
   else
       krange=$krange'33,55,71'
   fi
else
  krange=$(parameter krange)
fi
echo "krange for SPAdes:" $krange

spadespath=''
##default is 3.6.2
if [[ $(parameter spadesversion) == '3.8.1' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/tools/miniconda2/envs/spades3.8.1/bin
 echo 'SPAdes version: '$(parameter spadesversion)
elif [[ $(parameter spadesversion) == '3.6.2' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/bin
 echo 'SPAdes version: '$(parameter spadesversion)
elif [[ $(parameter spadesversion) == '3.5.0' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.5.0
 echo 'SPAdes version: '$(parameter spadesversion)
elif [[ $(parameter spadesversion) == '3.7.0' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.7.0-Linux/bin
 echo 'SPAdes version: '$(parameter spadesversion)
elif [[ $(parameter spadesversion) == '3.8.0' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.8.0-Linux/bin
 echo 'SPAdes version: '$(parameter spadesversion)
elif [[ $(parameter spadesversion) == '3.8.2' ]]; then
 spadespath=/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.8.2-Linux/bin
 echo 'SPAdes version: '$(parameter spadesversion)
else
 echo 'SPAdes version not found. Uses version 3.6.2'
 spadespath=/hpc/local/CentOS7/dla_mm/bin
fi

#echo ""
#echo "run spades"
#if [[ ! -f $assemblyfolder/scaffolds.fasta ]]; then
#$spadespath/spades.py $(parameter spadesparam) -1 "$DIR"/"$samplename".R1.fastq -2 "$DIR"/"$samplename".R2.fastq -k "$krange" --cov-cutoff $(parameter mincov) -o "$assemblyfolder"
#fi

#echo ""
#echo "run unicycler"
source activate /hpc/local/CentOS7/dla_mm/tools/miniconda2/envs/unicycler
unicycler -1 "$DIR"/"$samplename".R1.fastq -2 "$DIR"/"$samplename".R2.fastq -o  "$assemblyfolder"


#echo ""
#echo "rerun spades if not completed"
#while [[ ! -f $assemblyfolder/scaffolds.fasta ]]; do
 #   $spadespath/spades.py $(parameter spadesparam) -1 "$DIR"/"$samplename".R1.fastq -2 "$DIR"/"$samplename".R2.fastq -k "$krange" --continue --cov-cutoff $(parameter mincov) -o "$assemblyfolder"
#done

#fi

if [[ ! -f $assemblyfolder/assembly.fasta ]]; then
  echo $samplename "has no scaffolds associated" ; exit 0
fi

echo "end of assembly"
echo "#############################################################"
echo "post-processing starts here"

mkdir -p $resultsdir/scaffolds

RENAMED=$resultsdir/scaffolds/"$samplename".fna

if [[ ! -f $RENAMED ]]; then
echo "filtering of contigs shorter than minlen"
seqtk seq -L $(parameter minlen) $assemblyfolder/assembly.fasta | sed s/NODE/$samplename"_"$V/g > $RENAMED
fi

sleep 20


if [[ ! -f $RENAMED ]]; then
 echo $samplename "has no scaffolds associated" ; exit 0
fi


PROKKAPREFIX=$(echo $samplename| cut -d'-' -f 4)

if [[ $(parameter annotation) == yes ]]; then
 echo "run prokka"
    if [[ -z $genus ]];then
      prokka --outdir $assemblyfolder/annotation  --prefix $samplename --centre $PROKKAPREFIX $(parameter prokkaparam) $RENAMED
    else 
      prokka --outdir $assemblyfolder/annotation  --prefix $samplename --centre $PROKKAPREFIX $(parameter prokkaparam) --genus $genus --usegenus $RENAMED
    fi
fi

mkdir -p  $resultsdir/annotation
cp $assemblyfolder/annotation/*faa $resultsdir/annotation/ 
cp $assemblyfolder/annotation/*gbk $resultsdir/annotation/
cp $assemblyfolder/annotation/*gff $resultsdir/annotation/



mkdir -p $resultsdir/taxonomy
## run kraken
if [[ $(parameter taxonomy) == yes ]]; then
  echo "run taxonomy determination"
  krakenoutput=$assemblyfolder/"$samplename".krakenout
  if [[ ! -f $krakenoutput ]]; then
     kraken $(parameter krakenparam) \
      --output $krakenoutput --fasta-input $RENAMED
     cut -f2,3 $krakenoutput > $assemblyfolder/"$samplename".kronain
     ktImportTaxonomy $assemblyfolder/$samplename.kronain -o $resultsdir/taxonomy/"$samplename"_taxonomy.html
  fi
fi



echo "Quality check"
## genome coverage
#function to calculate genome coverage
genomesize () {
  grep -v ">" $1 | wc -m
}


if [[ $(parameter summary) == yes ]]; then
    qualitytext=$resultsdir/Coverage_"$project".tsv
    if [ -f $qualitytext ]; then
     cat $DIR/$samplename.*.fastq > $DIR/tmp$samplename.fastq
     seqtk fqchk $DIR/tmp$samplename.fastq | grep ALL | sed s/ALL/$samplename/g |\
     sed s/".fastq"//g |  awk '{ printf  "%s\t%.f\t%.f\n", $1, $2, $2/"'$(genomesize $RENAMED)'" }'>> $qualitytext
     rm $DIR/tmp$samplename.fastq
    else 
     echo -e "Sample\tNumber_of_bases\tEstimated_Coverage(x)" > $qualitytext
     cat $DIR/$samplename.*.fastq > $DIR/tmp$samplename.fastq
     seqtk fqchk $DIR/tmp$samplename.fastq | grep ALL | sed s/ALL/$samplename/g |\
     sed s/".fastq"//g |  awk '{ printf  "%s\t%.f\t%.f\n", $1, $2, $2/"'$(genomesize $RENAMED)'" }'>> $qualitytext
     rm $DIR/tmp$samplename.fastq
    fi
fi


