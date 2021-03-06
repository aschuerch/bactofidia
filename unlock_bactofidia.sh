#!/bin/bash -i
# command to unlock the workflow

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo 
echo "snakemake will be installed" 
conda install -y snakemake=5.14.0
fi


#remove all generated temporary directories if any
for i in tmp data stats results log
 do
 rm -r "$i"
 done

snakemake --snakefile Snakefile.assembly  --cores all --unlock --config configfile=config/config.yaml
