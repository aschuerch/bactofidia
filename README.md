### bactofidia: Basic microbial WGS analysis pipeline

*bactofidia* is a bacterial assembly and basic analysis pipeline using Snakemake and bioconda.

## Usage (short version)

Clone this repository with

```bash
git clone https://github.com/aschuerch/bactofidia.git
```

Move your paired-end read sequencing files (Sample1_R1.fastq.gz, Sample1_R2.fastq.gz, Sample1_R1.fastq.gz and Sample1_R2.fastq.gz) 
to this folder (or symlink them). Run the pipeline with


```bash
./bactofidia.sh Sample1 Sample2
```


## Prerequisites

Bactofidia runs under bash and relies heavily on a miniconda installation is required. If miniconda is not present, 
the script will attempt to install the latest version in a temporary folder.
However, I would recommend to install miniconda. See the [bioconda page](https://bioconda.github.io)

# Usage (long version)

The pipeline takes Illumina paired-end sequencing reads as compressed sequencing files (.fastq.gz) 
which must be present in the same folder from where the script is run.

Use only the sample name to call the script (not the .fastq.gz ending)

The config.yaml file can be adjusted for parameters of the different tools.

The different versions of the packages that are run are defined in the envs folder. 

For debugging or testing purposes, the pipeline itself can be dry-run with 

```bash
source activate snakemake
snakemake -np  --snakefile Snakefile.assembly --config configfile=config.yaml
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

Command line parameters for the different tools can be adjusted in the config.yaml or config_miseq.yaml file or 
in the Snakefile.assembly directly. For many cases, the default parameters should be sufficient


## Using different package versions

Package versions can be adjusted in envs/*yaml. 
The packages are in different files mainly due to different dependencies.
Please visit [bioconda](http://bioconda.github.io/) for available packages


## Adding other tools

For further customizing, see [snakemake documentation](https://snakemake.readthedocs.io/en/stable/)


## Trouble shooting

Sometimes the snakemake workflow will give you an error due to an already running instance.
In this case, unlock the snakemake instance with

```bash
./unlock.sh
```


