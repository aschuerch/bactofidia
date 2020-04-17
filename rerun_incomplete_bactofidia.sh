#!/bin/bash -i

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo 
echo "snakemake will be installed" 
conda install -y snakemake=5.14.0
fi

snakemake --snakefile Snakefile.assembly  --cores all --use-conda --printshellcmds  --config configfile=config/config.yaml --rerun-incomplete


