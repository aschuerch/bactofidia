#!/bin/bash -i

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo 
echo "snakemake will be installed" 
conda install -y snakemake=5.14.0
fi

# run snakemake just to unlock the folder

snakemake --snakefile Snakefile.assembly --cores all --use-conda  --printshellcmds --config configfile=config/config.yaml --unlock


## reuse log file

log=$(pwd)/log/call_assembly.txt


## re-run the snakemake pipeline

snakemake --snakefile Snakefile.assembly  --cores all --use-conda --printshellcmds  --config configfile=config/config.yaml --rerun-incomplete 2>&1 | tee -a "$log"
