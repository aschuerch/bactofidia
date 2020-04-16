#!/bin/bash -i

# Check if snakemake is found or install directly into base
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo
echo "snakemake will be installed"
conda install -y snakemake
fi

#dryrun
conda activate snakemake

snakemake -np  --snakefile Snakefile.assembly --config configfile=config/config.yaml

conda deactivate
