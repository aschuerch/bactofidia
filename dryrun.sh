#!/bin/bash
source activate snakemake
snakemake -np  --snakefile Snakefile.assembly --config configfile=config.yaml
source deactivate

