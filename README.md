### bactofidia: Basic microbial WGS analysis pipeline

*bactofidia* is a bacterial assembly and basic analysis pipeline using Snakemake and bioconda.

## Usage (short version)

Clone this repository with

```bash
git clone https://github.com/aschuerch/bactofidia.git
```

Move your paired-end read sequencing files (Sample1_R1.fastq.gz, Sample1_R2.fastq.gz, Sample1_R1.fastq.gz and Sample1_R2.fastq.gz) to this folder (or symlink them). Run the pipeline with


```bash
./bactofidia.sh Sample1 Sample2
```


## Installation

A bioconda installation is required.

This is an example script to install bioconda on a Linux 64-bit architecture. For other options, see the [bioconda page](https://bioconda.github.io)

```bash
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
```
Consider to add the miniconda directory to your PATH environment variable by editing your ~/.bash_profile

Next, create a virtual environment for snakemake using python3:

```bash
 conda create -y -n snakemake python=3.5 snakemake
```


## Usage (long version)

Before running the pipeline for the first time, another virtual environment needs to be created. Packages and versions are specified in 'package-list.txt'. Adjust this file to your needs. See bioconda.github.io for available packages.

Create the environment with 

```bash
conda create --file package-list.txt -n [bactofidia_custom]
```
where [bactiofidia_custom] matches the name of the environment given in the config.yaml

The config.yaml file can be adjusted for parameters of the different tools.

For debugging or testing purposes, the pipeline itself can be dry-run with 

```bash
source activate snakemake
snakemake -np --config configfile=config.yaml
source deactivate
```

The pipeline takes compressed sequencing files (.fastq.gz) which must be present in the same folder from where the script is called.
Use only the sample name to call the script (not the .fastq.gz ending).

## Testing

Test the whole pipeline with:

```bash
cp test/Test*gz .
./bactofidia.sh Test
```

This will run the pipeline on the included Test.fastq.gz files.

## Components

Currently it runs:
 - quality check before trimming using fastqc
 - trimming with [seqtk](http://bioconda.github.io/recipes/seqtk/README.html)
 - assembly with [spades](http://bioconda.github.io/recipes/spades/README.html)
 - mlst with [mlst](http://bioconda.github.io/recipes/mlst/README.html)
 - resistance gene determination with [abricate](http://bioconda.github.io/recipes/abricate/README.html) using the [resfinder database](https://cge.cbs.dtu.dk/services/ResFinder/)
 - quality assessment of assembly with [quast](http://bioconda.github.io/recipes/quast/README.html)
 - coverage estimation with [bbmap2](http://bioconda.github.io/recipes/bbmap/README.html)
 - summarizing report with [multiqc](http://bioconda.github.io/recipes/multiqc/README.html)

Running only 

```bash
./bactofidia.sh
```

will give an explanation of the (limited) options.


## Adjusting command line parameters

Command line parameters for the different tools can be adjusted in the config.yaml or config_miseq.yaml file or in Snakefile directly. For many cases, the default parameters should be sufficient


## Using different package versions

Package versions can be adjusted in package-list.txt. Please visit [bioconda](http://bioconda.github.io/) for available packages


## Adding other tools

For further customizing, see [snakemake documentation](https://snakemake.readthedocs.io/en/stable/)


## Trouble shooting

Unlock the snakemake instance
snakemake --snakefile Snakefile.assembly --unlock --config configfile=config.yaml
