import os
from snakemake.utils import min_version

min_version("3.9")
configfile: "config.yaml"

##to build a virt env
seqtkversion = config["seqtk"]["version"] 
spadesversion = config["SPAdes"]["version"]
versiontag = config["SPAdes"]["versiontag"]

# Prokka
prokkaversion = config["prokka"]["version"]
# mlst
mlstversion = config["mlst"]["version"]

# QUAST  
quastversion = config["QUAST"]["version"]

# Kraken
krakenversion = config["kraken"]["version"] 

# Krona
kronaversion = config["krona"]["version"]
# Checkm
checkmversion = config["checkm"]["version"]


# python 2 virtual environments
os.system ( "conda create -y -n bactofidia  python=2.7  checkm-genome={} quast={} seqtk={} spades={} prokka={} mlst={} kraken={} krona={}".format ( checkmversion, quastversion, seqtkversion, spadesversion, prokkaversion, mlstversion, krakenversion, kronaversion)    )

os.system ( "conda-env export -n bactofidia > bactofidia.yml") 


# Collect samples

SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)


onsuccess:
    # delete virtual environment
    os.system ( "conda-env remove -n bactofidia") 
    print("Workflow finished!")


onerror:
    # delete virtual environment
    os.system ( "conda-env remove -n bactofidia") 
    print("Workflow finished")



rule all:
    input:
       "stats/Trimmingstats.tsv",
       "stats/MLST.tsv",
       "stats/AssemblyQC.html",
       "stats/AssemblyQC.tsv",
       "stats/ResFinder.tsv",
       "stats/SpeciesDetermination.tsv",
       expand ("scaffolds/{sample}.fna", sample=SAMPLES),
       expand ("annotation/{sample}.gff", sample=SAMPLES),
       expand ("stats/Taxonomy_{sample}.html", sample=SAMPLES)


rule fastqc_before:
    input:
        "data/{sample}_{r}.fastq.gz"
    output:
        temp("stats/{sample}_{r}_Trimmingstats_before_trimming")
    conda:
        "bactofidia.yml"
    shell:
        "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"

rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("trimmed/{sample}_{r}.fastq")
     params:
        p = config["seqtk"]["params"]
     conda:
        "bactofidia.yml"
     shell: 
        "seqtk trimfq {params.p} {input} > {output}"

rule fastqc_after:
    input:
        "trimmed/{sample}_{r}.fastq"
    output:
         temp("stats/{sample}_{r}_Trimmingstats_after_trimming")
    conda:
        "bactofidia.yml"
    shell:
        "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"


rule trimstat:
    input:
        after=expand("stats/{sample}_{r}_Trimmingstats_after_trimming", sample=SAMPLES, r=R),
        before=expand("stats/{sample}_{r}_Trimmingstats_before_trimming", sample=SAMPLES, r=R)
    output:
        "stats/Trimmingstats.tsv"
    conda:
        "bactofidia.yml"
    shell:
        "echo -e 'Reads\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high' > {output} ;"
        "cat {input.before} >> {output} ;"
        "cat {input.after} >> {output} ;"
       
rule spades:
    input: 
        R1="trimmed/{sample}_R1.fastq",
        R2="trimmed/{sample}_R2.fastq"
    output:
        "assembly/{sample}/scaffolds.fasta"
    params:
        spadesparams = config["SPAdes"]["params"],
        kmer = config.get("krange"),  #if else statement desired
        #kmer = config["krange"],
        cov = config["mincov"],
       # spadesversion = spadesversion,
        outfolder = "assembly/{sample}"
    conda:
        "bactofidia.yml"
    shell:
        "spades.py -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"

rule rename:
    input:
        "assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna"
    params:
        minlen = config["minlen"],
        versiontag = "{sample}:"+versiontag
    conda:
        "bactofidia.yml"
    shell:
        "seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.versiontag}/g > {output}"

rule annotation:
    input:
        "scaffolds/{sample}.fna"
    output:
        "annotation/{sample}.gff"
    params:
        dir = "annotation",
        params = config["prokka"]["params"],
        prefix = "{sample}"
    conda:
        "bactofidia.yml"
    shell:
        "prokka --force --prefix {params.prefix} --outdir {params.dir} {params.params} {input} "

        
rule taxonomy_1:
    input:
        "scaffolds/{sample}.fna"
    output:
        temp("taxonomy/{sample}.krakenout")
    params:
        p=config["kraken"]["params"]
    conda:
        "bactofidia.yml"
    shell:
        "kraken {params.p} --output {output} --fasta_input {input}"
        "source deactivate"

rule taxonomy_2:
    input:
        "taxonomy/{sample}.krakenout"
    output:
        temp("taxonomy/{sample}.kronain")
    conda:
        "bactofidia.yml"
    shell:
        "cut -f2,3 {input} > {output}"


rule taxonomy_3:
    input:
        "taxonomy/{sample}.kronain"
    output:
        "stats/Taxonomy_{sample}.html"
    conda:
        "bactofidia.yml"
    shell:
        "ktImportTaxonomy {input} -o {output}"

rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    params:
        config["mlst"]["params"]
    conda:
        "bactofidia.yml"
    shell:
        "mlst {params} {input} >> {output}"


rule quast:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        html = "stats/AssemblyQC.html",
        tsv = "stats/AssemblyQC.tsv"
    params:
        outfolder = "stats/quasttemp",
        p = config["QUAST"]["params"]
    conda:
        "bactofidia.yml"
    shell:
        "quast.py {params.p} -o {params.outfolder} {input}"
        " && mv stats/quasttemp/report.html {output.html}"
        " && mv stats/quasttemp/report.tsv {output.tsv}"
        " && rm -r stats/quasttemp"


rule resfinder:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES) 
    output:
        "stats/ResFinder.tsv"
    conda:
        "bactofidia.yml"
    shell:
        "abricate {input} > {output}" 
        "&& sed -i 's/scaffolds\///g' {output}"


rule checkm:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        file = "stats/SpeciesDetermination.tsv",
        folder = temp("checkm")
    params:
        config["checkm"]["params"]
    conda:
        "bactofidia.yml"
    shell:
        "checkm lineage_wf scaffolds {output.folder} {params} --tab_table -x fna > {output.file}"

