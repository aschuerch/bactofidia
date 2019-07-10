#!/bin/bash
# command to unlock the workflow
source activate snakemake

#remove all generated temporary directories if any
for i in tmp data stats results log
 do
 rm -r "$i"
 done


snakemake --snakefile Snakefile.assembly --unlock --config configfile=config/config.yaml
source deactivate
