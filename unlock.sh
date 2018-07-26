#!/bin/bash
# command to unlock the workflow
source activate snakemake
snakemake --snakefile Snakefile.assembly --unlock --config configfile=config.yaml
source deactivate
