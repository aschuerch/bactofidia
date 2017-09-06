import os
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
    # delete files 
    print("Workflow finished!")


onerror:
    print("Workflow finished with errors")


rule all:
    input:
        "stats/multiqc_report.html",
        "stats/AssemblyQC.html",
        "stats/ResFinder.tsv",
        "stats/MLST.tsv",
#        expand ("data/{sample}_R1.fastq.gz.msh", sample=SAMPLES),
         expand ("scaffolds/{sample}.fna", sample=SAMPLES),
        expand ("stats/CoverageStatistics_{sample}_scaffolds.txt",sample=SAMPLES),
        expand("stats/CoverageStatistics_{sample}.txt", sample=SAMPLES),
#        expand ("annotation/{sample}.gff", sample=SAMPLES),
#        expand ("stats/Taxonomy_{sample}.html", sample=SAMPLES),
#        expand ("stats/Trimmingstats_{sample}.tsv", sample=SAMPLES),
#       "stats/Trimmingstats.tsv"


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

rule multiqc:
    input:
        expand("tmp/{sample}_R1_fastqc.html", sample=SAMPLES)
    output:
        "stats/multiqc_report.html"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u "
        "&& multiqc tmp/ -n {output}"
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

#rule sketch:
 #   input:
  #      "data/{sample}_R1.fastq.gz"
   # output:
    #    "data/{sample}_R1.fastq.gz.msh" 
   ## params:
     #   p = config["mash"]["params"],
      #  virtenv = config["virtual_environment"]["name"]
   # shell:
    ##    "set +u; source activate {params.virtenv}; set -u"
      #  "&& mash sketch {params.p} {input}"
       # "&& set +u; source deactivate; set -u"
   

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
        
#rule mash_dist:
 #   input:
  #      query = "data/{sample}_R1.fastq.gz,
   #     ref = "data/ref.msh"
   # output:
    #    "ref/name_{sample}.txt"
   # params:
    ##    p = config["mash"]["params"],
      #  virtenv = config["virtual_environment"]["name"]
   # shell:
    #    "set +u; source activate {params.virtenv}; set -u"
     #   "&& mash dist -t {params.p} {input.query} {input.ref} | sort -k2 | head -n 1 | cut -f 1,1 | cut -f 1,1 -d / >> {output} "

#rule get_ref:
 #   input:
  #      "ref/name_{sample}.txt"
   # output:
    #    url = "ref/url_{sample}.txt",
     #  seq = "ref/refseq_{sample}.txt",
      #  dest = "ref/{sample}
    #shell:
     #   "echo -n 'wget -P {output.dest}' >> {output.url}"
      #  "&&grep  -f {input} ftpdirpaths | tr -d '\n' >> {output.url}"
       # "&&echo '/_genomic.fna.gz' >> {output.url}"
        #"&& bash {output.url}"


#rule taxonomy_1:
 #   input:
  #      "scaffolds/{sample}.fna"
   # output:
    #    temp("tmp/{sample}.krakenout")
    #params:
     #   p=config["kraken"]["params"],
      #  virtenv = config["virtual_environment"]["name"]
   # shell:
    #    "set +u; source activate {params.virtenv}; set -u"
     #   "&& kraken {params.p} --output {output} --fasta_input {input}"
      #  "&& set +u; source deactivate; set -u"

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
        html = "stats/AssemblyQC_{sample}.html"
    params:
        outfolder = "stats/quast",
        p = config["QUAST"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& quast {params.p} -o {params.outfolder} {input}"
        "&& cp stats/quast/report.html {output.html}"
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
        ref = expand("scaffolds/{sample}.fna", sample=SAMPLES),
        in_R1 = expand("tmp/{sample}_R1.fastq", sample=SAMPLES, r=R),
        in_R2 = expand("tmp/{sample}_R2.fastq", sample=SAMPLES, r=R)
    output:
        concat = temp("tmp/{sample}.fastq"),
        covstats_detail = "stats/CoverageStatistics_{sample}_scaffolds.txt",
        covstats = "stats/CoverageStatistics_{sample}.txt"
    params:
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& cat {input.in_R1} {input.in_R2} > {output.concat}"
        "&& bbmap.sh in={output.concat} ref={input.ref} covstats={output.covstats_detail} >> {output.covstats} 2>&1 " 
        "&& set +u; source deactivate; set -u"


#rule stat:
 #   input:
  #      R1before = "data/{sample}_R1.fastq.gz",
   #     R1after = "tmp/{sample}_R1.fastq", 
    #    R2before = "data/{sample}_R2.fastq.gz",
     #   R2after = "tmp/{sample}_R2.fastq", 
#    output:
 #       "stats/Trimmingstats_{sample}.tsv"
 
