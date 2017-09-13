import os
from snakemake.utils import min_version

#min_version("3.9")
configfile : config.get("configfile") 

def kmer_determination():
    if (config.get("krange")):
        kmer = config.get("krange")
    else:
        kmer = config["SPAdes"]["krange"]
    return kmer

versiontag = config["virtual_environment"]["name"]

print(versiontag)

# Collect samples
SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)
R = set(R)

onsuccess:
    # delete files 
    # os.remove (tmp)
    print("Workflow finished!")


onerror:
    print("Workflow finished with errors")


rule all:
    input:
        "stats/multiqc_report.html",
#        "stats/ResFinder.tsv",
#        "stats/MLST.tsv",
#        "stats/CoverageStatistics_summary.tsv",
        expand ("scaffolds/{sample}.fna", sample=SAMPLES)
#        expand ("stats/CoverageStatistics_{sample}_scaffolds.txt",sample=SAMPLES),
#        expand("stats/CoverageStatistics_{sample}.txt", sample=SAMPLES)


rule fastqc:
    input:
        "data/{sample}_R1.fastq.gz"
    output:
        temp("tmp/{sample}_R1_fastqc.html"),
        temp = temp("data/{sample}_R1.fastq")
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u "
        "&& zcat {input} >> {output.temp} "
        "&& fastqc {output.temp} -o tmp"
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
        outfolder = "tmp/assembly/{sample}",
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
        expand("scaffolds/{sample}.fna", sample=SAMPLES)      
    output:
        html = temp("tmp/report.html")
    params:
        outfolder = "tmp",
        p = config["QUAST"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& quast {params.p} -o {params.outfolder} {input}"
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

rule map:
    input:
        ref = "scaffolds/{sample}.fna",
        in_R1 = "tmp/{sample}_R1.fastq",
        in_R2 = "tmp/{sample}_R2.fastq"
    output:
        concat = temp("tmp/{sample}_concatenated.fastq"),
        covstats_detail = "stats/CoverageStatistics_{sample}_scaffolds.txt",
        covstats = "stats/CoverageStatistics_{sample}.txt"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& cat {input.in_R1} {input.in_R2} > {output.concat}"
        "&& bbmap.sh in={output.concat} ref={input.ref} covstats={output.covstats_detail} >> {output.covstats} 2>&1 " 
        "&& set +u; source deactivate; set -u"

rule sum:
    input:
        expand("stats/CoverageStatistics_{sample}.txt", sample=SAMPLES)
    output:
        "stats/CoverageStatistics_summary.tsv"
    shell:
        "grep 'Average coverage' {input} | sed 's@CoverageStatistics\_@@g' | sed 's@\:Average coverage\:@@g' >> {output}"


rule multiqc:
    input:
        expand("tmp/{sample}_R1_fastqc.html", sample=SAMPLES),
        "tmp/report.html"
    output:
        "stats/multiqc_report.html"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u "
        "&& multiqc tmp/ -n {output}"
        "&& set +u; source deactivate; set -u"

