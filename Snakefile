from snakemake.utils import min_version
min_version("3.2")
configfile: "config.yaml"



seqtk_version = config["seqtk"]["version"] 
seqtkparams = config["seqtk"]["params"]

#Assembly with SPAdes
spadesparams = config["SPAdes"]["params"]
spadesversion = config["SPAdes"]["version"]
krange = config.get("krange")
#krange = "57,97,127" #adjust if you would like to select a different kmer range for spades assembly
versiontag = config["SPAdes"]["versiontag"]

#Minimum length of contigs
minlen = config["minlen"]
#Minimum coverage depth of contigs
mincov = config["mincov"]

# QUAST  
quastparams = config["QUAST"]["params"]

# Prokka
prokkaparams = config["prokka"]["params"] 

# Kraken
krakenparams = config["kraken"]["params"] 

# Checkm
checkmparams = config["checkm"]["params"] 


mlstparams = config["mlst"]["params"]


SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)

def determine_spadespath(spadesversion): #some are not available as conda package
    spadespath = "/hpc/local/CentOS7/dla_mm/bin/spades.py" #standard
    if spadesversion == "3.9.0":
        spadespath = "/hpc/local/CentOS7/dla_mm/tools/miniconda2/envs/snakemake-tutorial/bin/spades.py"
    if spadesversion == "3.8.2":
        spadespath = "/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.8.2-Linux/bin/spades.py"
    if spadesversion == "3.8.1":
        spadespath = "/hpc/local/CentOS7/dla_mm/tools/miniconda2/envs/spades3.8.1/bin/spades.py"
    if spadesversion == "3.8.0":  
        spadespath = "/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.8.0-Linux/bin/spades.py"
    if spadesversion == "3.6.2":
        spadespath = "/hpc/local/CentOS7/dla_mm/bin/spades.py"
    if spadesversion == "3.5.0":
        spadespath = "/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.5.0/spades.py"
    if spadesversion == "3.7.0":
        spadespath = "/hpc/local/CentOS7/dla_mm/bin/SPAdes-3.7.0-Linux/bin/spades.py"
    return spadespath

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
    shell:
        "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"

rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("trimmed/{sample}_{r}.fastq")
     params:
        seqtkparams
     shell: 
        "seqtk trimfq {params} {input} > {output}"


rule fastqc_after:
    input:
        "trimmed/{sample}_{r}.fastq"
    output:
         temp("stats/{sample}_{r}_Trimmingstats_after_trimming")
    shell:
         "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"

rule trimstat:
    input:
        after=expand("stats/{sample}_{r}_Trimmingstats_after_trimming", sample=SAMPLES, r=R),
        before=expand("stats/{sample}_{r}_Trimmingstats_before_trimming", sample=SAMPLES, r=R)
    output:
        "stats/Trimmingstats.tsv"
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
        spadesparams = spadesparams,
        kmer = krange,
        cov = mincov,
        spadesversion = determine_spadespath(spadesversion),
        outfolder = "assembly/{sample}"
    shell:
        "python {params.spadesversion} -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"

rule rename:
    input:
        "assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna"
    params:
        minlen = minlen,
        versiontag = "{sample}:"+versiontag
    shell:
        "seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.versiontag}/g > {output}"

rule annotation:
   input:
        "scaffolds/{sample}.fna"
   output:
        "annotation/{sample}.gff"
   params:
        dir = "annotation",
        params = prokkaparams,
        prefix = "{sample}"
   shell:
        "prokka --force --prefix {params.prefix} --outdir {params.dir} {params.params} {input} "
        
rule taxonomy_1:
   input:
        "scaffolds/{sample}.fna"
   output:
        temp("taxonomy/{sample}.krakenout")
   params:
        krakenparams
   shell:
        "kraken {params} --output {output} --fasta_input {input}"

rule taxonomy_2:
   input:
        "taxonomy/{sample}.krakenout"
   output:
        temp("taxonomy/{sample}.kronain")
   shell:
        "cut -f2,3 {input} > {output}"

rule taxonomy_3:
   input:
        "taxonomy/{sample}.kronain"
   output:
        "stats/Taxonomy_{sample}.html"
   shell:
        "ktImportTaxonomy {input} -o {output}"


rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    params:
        mlstparams
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
        minlen = minlen
    shell:
        "quast.py --min-contig {params.minlen} -o {params.outfolder} {input}"
        " && mv stats/quasttemp/report.html {output.html}"
        " && mv stats/quasttemp/report.tsv {output.tsv}"
        " && rm -r stats/quasttemp"
 
rule resfinder:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES) 
    output:
        "stats/ResFinder.tsv"
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
       checkmparams
    shell:
       "set +u; source activate checkm; set -u;"
       "checkm lineage_wf scaffolds {output.folder} {params} --tab_table -x fna > {output.file} ;"
       "source deactivate"
