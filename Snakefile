import os
import sys
import glob
import subprocess
import shutil
from random import randint


from snakemake.utils import min_version

min_version("3.9")
configfile: "config.yaml"


def kmer_determination():
    if (config.get("krange")):
        kmer = config.get("krange")
    else:
        kmer = config["SPAdes"] ["krange"]
    return kmer

versiontag = config["virtual_environment"]["name"]

print(versiontag)

# Collect samples
SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)
R = set(R)

onsuccess:
    # delete virtual environment
#    for i in virtenvs:
 #       os.system ( "conda-env remove -y -n {}".format (i)) 
    print("Workflow finished!")


onerror:
 #   # delete virtual environment
  #  for i in virtenvs:
   #     os.system ( "conda-env remove -y -n {}".format(i)) 
    print("Workflow finished")


rule all:
    input:
       "stats/Trimming.tsv",
       "stats/MLST.tsv",
       "stats/AssemblyQC.html",
       "stats/AssemblyQC.tsv",
       "stats/ResFinder.tsv",
       "stats/MLST.tsv",
       expand ("scaffolds/{sample}.fna", sample=SAMPLES),
       expand ("annotation/{sample}.gff", sample=SAMPLES),
       expand ("stats/Taxonomy_{sample}.html", sample=SAMPLES),


rule fastqc_before:
    input:
        "data/{sample}_{r}.fastq.gz"
    output:
        temp("tmp/{sample}_{r}_Trimmingstats_before_trimming")
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& echo -ne {input}'\t'| sed -n '2~4p' {input} | wc -m >> {output}"
        "&& set +u; source deactivate; set -u"

rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("tmp/{sample}_{r}.fastq")
     params:
        p = config["seqtk"]["params"],
        virtenv = config["virtual_environment"]["name"]
     shell: 
        "set +u; source activate {params.virtenv}; set -u"
        "&& seqtk trimfq {params.p} {input} > {output} "
        "&& set +u; source deactivate; set -u"

rule fastqc_after:
    input:
        temp("tmp/{sample}_{r}.fastq")
    output:
        temp("tmp/{sample}_{r}_Trimmingstats_after_trimming")
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& echo -ne {input}'\t'| sed -n '2~4p' {input} | wc -m >> {output}"
        "&& set +u; source deactivate; set -u"       

rule trimstat:
    input:
        after=expand("tmp/{sample}_{r}_Trimmingstats_after_trimming", sample=SAMPLES, r=R),
        before=expand("tmp/{sample}_{r}_Trimmingstats_before_trimming", sample=SAMPLES, r=R)
    output:
        "stats/Trimming.tsv"
    shell:
        "echo -e 'Reads\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high' > {output} "
        "&& echo -en {input.before}  >> {output} "
        "&& cat {input.before} >> {output} "
        "&& echo -en {input.after} >> {output} "
        "&& cat {input.after} >> {output} "
       
rule spades:
    input: 
        R1=temp("tmp/{sample}_R1.fastq"),
        R2=temp("tmp/{sample}_R2.fastq")
    output:
        temp("tmp/assembly/{sample}/scaffolds.fasta")
    params:
        spadesparams = config["SPAdes"]["params"],
        kmer = kmer_determination(),
        cov = config["mincov"],
        outfolder = "assembly/{sample}",
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& spades.py -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"
        "&& set +u; source deactivate; set -u"

rule rename:
    input:
        temp("tmp/assembly/{sample}/scaffolds.fasta")
    output:
        "scaffolds/{sample}.fna",
    params:
        minlen = config["minlen"],
        versiontag = "{sample}_"+versiontag,
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.versiontag}/g > {output}"
        "&& set +u; source deactivate; set -u"
      
rule annotation:
    input:
        "scaffolds/{sample}.fna"
    output:
        "annotation/{sample}.gff",
    params:
        dir = "annotation",
        params = config["prokka"]["params"],
        prefix = "{sample}",
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& prokka --force --prefix {params.prefix} --outdir {params.dir} {params.params} {input} "
        "&& set +u; source deactivate; set -u"
        
rule taxonomy_1:
    input:
        "scaffolds/{sample}.fna"
    output:
        temp("tmp/{sample}.krakenout")
    params:
        p=config["kraken"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& kraken {params.p} --output {output} --fasta_input {input}"
        "&& set +u; source deactivate; set -u"

rule taxonomy_2:
    input:
        "tmp/{sample}.krakenout"
    output:
        temp("tmp/{sample}.kronain")
    shell:
        "cut -f2,3 {input} > {output}"

rule taxonomy_3:
    input:
        temp("tmp/{sample}.kronain")
    output:
        "stats/Taxonomy_{sample}.html"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& ktImportTaxonomy {input} -o {output}"
        "&& set +u; source deactivate; set -u"

rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    params:
        mlst = config["mlst"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& mlst {params.mlst} {input} > {output}"
        "&& set +u; source deactivate; set -u"

rule quast:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES),
        
    output:
        html = "stats/AssemblyQC.html",
        tsv = "stats/AssemblyQC.tsv"
    params:
        outfolder = "stats/quasttemp",
        p = config["QUAST"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& quast {params.p} -o {params.outfolder} {input}"
        "&& mv stats/quasttemp/report.html {output.html}"
        "&& mv stats/quasttemp/report.tsv {output.tsv}"
        "&& rm -r stats/quasttemp"
        "&& set +u; source deactivate; set -u"

rule resfinder:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES) 
    output:
        "stats/ResFinder.tsv"
    params:
        p = config["abricate"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& abricate {params.p} {input} > {output}" 
        "&& sed -i 's/scaffolds\///g' {output}"
        "&& set +u; source deactivate; set -u"

