import shutil

from snakemake.utils import min_version

min_version("3.9")

configfile : config.get("configfile") 

# Collect samples. They need to be stored in data/
SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)
R = set(R)

onsuccess:
    # delete files
    if  (config["remove_temp"] == "yes"):
        shutil.rmtree ("tmp")
    print("Workflow finished!")


onerror:
    print("Workflow finished with errors")


rule all:
    input:
        "stats/multiqc_report.html",
        "stats/ResFinder.tsv",
        "stats/MLST.tsv",
        "stats/CoverageStatistics_summary.tsv",
        expand ("scaffolds/{sample}.fna", sample=SAMPLES)


rule fastqc:
    input:
        expand("data/{sample}_R1.fastq.gz", r=R)
    output:
        "tmp/{sample}_R1_fastqc.html",
        temp = "data/{sample}_R1.fastq"
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
        "tmp/{sample}_{r}.fastq"
    params:
        p = config["seqtk"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell: 
        "set +u; source activate {params.virtenv}; set -u"
        "&& seqtk trimfq {params.p} {input} > {output} "
        "&& set +u; source deactivate; set -u"   

rule spades:
    input: 
        R1="tmp/{sample}_R1.fastq",
        R2="tmp/{sample}_R2.fastq"
    output:
        "tmp/assembly/{sample}/scaffolds.fasta"
    params:
        spadesparams = config["SPAdes"]["params"],
        kmer = config["SPAdes"]["krange"],
        cov = config["mincov"],
        outfolder = "tmp/assembly/{sample}",
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& spades.py -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"
        "&& set +u; source deactivate; set -u"

rule rename:
    input:
        "tmp/assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna",
    params:
        minlen = config["minlen"],
        versiontag = "{sample}_"+config["virtual_environment"]["name"],
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
        html = "tmp/report.html",
        html2 = "stats/Assembly_report.html"
    params:
        outfolder = "tmp/",
        p = config["QUAST"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& quast {params.p} -o {params.outfolder} {input}"
        "&& cp {output.html} {output.html2} " 
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
        covstats_detail = "tmp/CoverageStatistics_{sample}_scaffolds.txt",
        covstats = "stats/CoverageStatistics_{sample}.txt"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& bbmap.sh in={input.in_R1} in2={input.in_R2} ref={input.ref} covstats={output.covstats_detail} >> {output.covstats} 2>&1 " 
        "&& set +u; source deactivate; set -u"

rule sum:
    input:
        expand("stats/CoverageStatistics_{sample}.txt", sample=SAMPLES)
    output:
        "stats/CoverageStatistics_summary.tsv"
    shell:
        "grep 'Average coverage' {input} | sed 's@CoverageStatistics\_@@g' | sed 's@\:Average coverage\:@@g' | sed 's@stats@@g' | sed 's@\\\@@g'| sed 's@\.txt@@g' >> {output}"


rule multiqc:
    input:
        expand("tmp/{sample}_R1_fastqc.html", sample=SAMPLES),
        "tmp/report.html",
        expand("stats/CoverageStatistics_{sample}.txt", sample=SAMPLES)
    output:
        "stats/multiqc_report.html",
        temp
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u "
        "&& multiqc tmp/ -n {output}"
        "&& set +u; source deactivate; set -u"

