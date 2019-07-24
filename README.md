# bactofidia
                                                                                                                                                                             

<img src="noun_Snake_1753915_B30083.png" width="300px" style="display: block; margin: auto;" />


### Basic microbial WGS analysis pipeline

*bactofidia* is a bacterial assembly and basic analysis pipeline using Snakemake and bioconda.
It is currently written for paired-end Illumina data with length 250 or 150. The pipeline is written to ensure reproducibility, 
and creates virtual software environments with the software versions that are used for analysis.

## Dependencies

bactofidia runs under bash and relies on software available from [bioconda](https://bioconda.github.io) and a (mini)conda (<4.7) installation. 
If conda is not present, the script will attempt to install a compatible miniconda version in a temporary directory.

## Usage 

Clone this repository with

```bash
git clone https://gitlab.com/aschuerch/bactofidia.git bactofidia_[myproject]
```

where [myproject] is the name of your project.

Copy or symlink your paired-end read sequencing files
(Sample1_R1.fastq.gz, Sample1_R2.fastq.gz, Sample2_R1.fastq.gz and Sample2_R2.fastq.gz) to the bactofidia_[myproject] directory.
After succesful execution of the pipeline, these files will be removed from this folder. Make sure they have been stored elsewhere.


The first underscore in the sample names is regarded as the delimiter for the sample name.
Avoid other underscores in the samplenames.

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

The first time bactofidia is run, it generates all virtual environments which can take a considerable time depending on the speed of your internet connection.
Do not interrupt this process!


## De-bugging and testing

For debugging or testing purposes, the pipeline itself can be dry-run with 

```bash
./dryrun_bactofidia.sh
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
results/
|
├── config
│   ├── ... # contains all configuration definitions, e.g. parameter choices
│ 
├── envs
│   ├── ... # contains all environment definitions, e.g. used versions of programs
│ 
├── scaffolds
│   ├── Test.fna   # scaffolds
│   
└── stats
    ├── annotated  # results of PROKKA
    │   ├── Test
    │      ├── Test.err
    │      ├── Test.faa
    │      ├── Test.ffn
    │      ├── Test.fna
    │      ├── Test.fsa
    │      ├── Test.gbk
    │      ├── Test.gff
    │      ├── Test.log
    │      ├── Test.sqn
    │      ├── Test.tbl
    │      ├── Test.tsv
    │      └── Test.txt
    │   
    ├── CoverageStatistics_summary.tsv # Summary of coverage for all sample
    ├── Extra
    │   ├── CoverageStatistics_Test.txt # Coverage details on each sample
    │ 
    ├── MLST.tsv # MLST results in table format
    ├── MultiQC_report_data.zip # Contains data to generate the MultiQC report
    ├── MultiQC_report.html # Check this file first for overview of quality
    └── ResFinder.tsv # Resistance genes per sample in table format
```

Opening the MultiQC_report.html in a browser is a good starting point to judge quality of data and assembly.


## Adjusting command line parameters

Command line parameters for the different tools can be adjusted in the config.yaml or config_miseq.yaml file in the config/ directory, 
or in the Snakefile.assembly directly. For many cases, the default parameters should be sufficient.


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
./unlock_bactofidia.sh
```
The pipeline can be tested with

```bash
./dryrun_bactofidia.sh
```
