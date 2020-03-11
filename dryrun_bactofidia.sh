#!/bin/bash -i

#activate snakemake or generate a snakemake env if not present
conda activate snakemake \
|| conda create -y -n snakemake snakemake python=3.5 \
&& conda activate snakemake

#dryrun
snakemake -np  --snakefile Snakefile.assembly --config configfile=config/config.yaml

conda deactivate
