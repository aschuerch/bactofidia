#!/bin/bash -i

# Check if snakemake is found or install directly into base 
if command -v snakemake > /dev/null; then ##version?
echo "snakemake found"
else
echo 
echo "snakemake will be installed" 
conda install -y snakemake
fi

#if command -v qstat > /dev/null; then

#snakemake \
#  --snakefile Snakefile.assembly \
#  --latency-wait 60 \
#  --config configfile=config/config.yaml \
#  --rerun-incomplete \
#  --verbose \
#  --forceall \
#  --keep-going \
#  --restart-times 5\
#  --use-conda \
#  --cluster \
#  'qsub -V -cwd -l h_vmem=32G -l h_rt=04:00:00 -e log/ -o log/ ' \
#  --jobs 100

#else
snakemake --snakefile Snakefile.assembly --use-conda --printshellcmds  --config configfile=config/config.yaml --rerun-incomplete

#fi

conda deactivate

