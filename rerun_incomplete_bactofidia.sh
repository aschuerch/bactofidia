#!/bin/bash -i

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo 
echo "snakemake will be installed" 
conda install -y snakemake
fi

conda activate snakemake

snakemake --snakefile Snakefile.assembly --use-conda --printshellcmds  --config configfile=config/config.yaml --rerun-incomplete

conda deactivate

