# bactofidia
                                                                                                                                                                             

<img src="noun_Snake_1753915_B30083.png" width="300px" style="display: block; margin: auto;" />


### Basic microbial WGS analysis pipeline

*bactofidia* is a bacterial assembly and basic analysis pipeline using Snakemake and bioconda.
It is currently written for paired-end Illumina data with length 250 or 150. The pipeline is written to ensure reproducibility, 
and creates virtual software environments with the software versions that are used for analysis.

## Dependencies

bactofidia runs under bash and relies on software available from [bioconda](https://bioconda.github.io) and a conda installation. 
If conda is not present, 
the script will attempt to install the latest miniconda version in a temporary directory.

## Usage 

Clone this repository with

```bash
git clone git@gitlab.com:aschuerch/bactofidia.git bactofidia_[myproject]
```

where [myproject] is the name of your project.

Move your paired-end read sequencing files (Sample1_R1.fastq.gz, Sample1_R2.fastq.gz, Sample2_R1.fastq.gz and Sample2_R2.fastq.gz) to this directory, or symlink them. 
The first underscore is regarded as the delimiter of the sample name. It's advisable to rename samples with underscores in their samplenames.
Run the pipeline with


```bash
./bactofidia.sh Sample1_R1.fastq.gz Sample1_R2.fastq.gz Sample2_R1.fastq.gz Sample2_R2.fastq.gz
```
or 

```bash
./bactofidia.sh ALL
```

The pipeline takes Illumina paired-end sequencing reads as compressed sequencing files (.fastq.gz) 
which must be present in the same directory from where the script is run.

The config.yaml or config_miseq.yaml files in the config/ directory can be adjusted for parameters of the different tools.

The different versions of the packages that are run are defined in the `envs/` directory. 

The first time bactofidia is run it generates all virtual environments which can take a considerable time depending on the speed of your internet connection.
Do not interrupt this process!


## De-bugging and testing

For debugging or testing purposes, the pipeline itself can be dry-run with 

```bash
./dryrun.sh
```

The pipeline takes compressed sequencing files (.fastq.gz) which must be 
present in the same directory from where the script is called.

Test the whole pipeline with:

```bash
ln -s  test/Test*gz .
./bactofidia.sh Test_R1.fastq.gz Test_R2.fastq.gz
```

This will run the pipeline on the included test files.

## Analysis steps

Currently it runs:
 - quality check before trimming using [fastqc](http://bioconda.github.io/recipes/fastqc/README.html)
 - trimming with [seqtk](http://bioconda.github.io/recipes/seqtk/README.html)
 - assembly with [spades](http://bioconda.github.io/recipes/spades/README.html)
 - mlst with [mlst](http://bioconda.github.io/recipes/mlst/README.html)
 - resistance gene determination with [abricate](http://bioconda.github.io/recipes/abricate/README.html) using the [resfinder database](https://cge.cbs.dtu.dk/services/ResFinder/)
 - annotation (general, not genus-specific) with [prokka](http://bioconda.github.io/recipes/prokka/README.html)
 - quality assessment of assembly with [quast](http://bioconda.github.io/recipes/quast/README.html)
 - coverage estimation with [bbmap2](http://bioconda.github.io/recipes/bbmap/README.html)
 - summarizing report with [multiqc](http://bioconda.github.io/recipes/multiqc/README.html)

Running only 

```bash
./bactofidia.sh
```

will give an explanation of the (limited) options.


## Output

The output can be found in the 'results' directory which contains the following files representing the output of the different tools

```
.
├── scaffolds
│   └── Test.fna
└── stats
    ├── CoverageStatistics_summary.tsv
    ├── Extra
    │   ├── Assembly_report.html
    │   └── CoverageStatistics_Test.txt
    ├── MLST.tsv
    ├── MultiQC_report_data.zip
    ├── MultiQC_report.html
    └── ResFinder.tsv
```

Opening the MultiQC_report.html in a browser is a good starting point to judge quality of data and assembly.


## Adjusting command line parameters

Command line parameters for the different tools can be adjusted in the config.yaml or config_miseq.yaml file or 
in the Snakefile.assembly directly. For many cases, the default parameters should be sufficient.


## Using different package versions

Package versions can be adjusted in envs/*yaml. 
The packages are in different files, mainly due to different dependencies such as python 2 / python 3.
Please visit [bioconda](http://bioconda.github.io/) for available packages.


## Adding other tools

For further customizing, see [snakemake documentation](https://snakemake.readthedocs.io/en/stable/)


## Trouble shooting

Sometimes the snakemake workflow will give you an error due to an already running instance.
In this case, unlock the snakemake instance with

```bash
./unlock.sh
```
The pipeline can be tested with

```bash
./dryrun.sh
```

