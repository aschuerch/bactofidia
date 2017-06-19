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
       "stats/AssemblyQC.html",
       "stats/ResFinder.tsv",
       expand ("scaffolds/{sample}.fna", sample=SAMPLES),
       expand ("annotation/{sample}.gff", sample=SAMPLES),
       expand ("stats/Taxonomy_{sample}.html", sample=SAMPLES),
       expand ("stats/Summary_{sample}.tsv", sample=SAMPLES),
       "stats/Summary.tsv"


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
        temp("tmp/MLST.tsv")
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
        txt = temp("tmp/AssemblyQC.txt")
    params:
        outfolder = "stats/quasttemp",
        p = config["QUAST"]["params"],
        virtenv = config["virtual_environment"]["name"]
    shell:
        "set +u; source activate {params.virtenv}; set -u"
        "&& quast {params.p} -o {params.outfolder} {input}"
        "&& mv stats/quasttemp/report.html {output.html}"
        "&& mv stats/quasttemp/transposed_report.txt {output.txt}"
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

rule stat:
    input:
        R1before = "tmp/{sample}_R1_Trimmingstats_before_trimming",
        R1after = "tmp/{sample}_R1_Trimmingstats_after_trimming", 
        R2before = "tmp/{sample}_R2_Trimmingstats_before_trimming",
        R2after = "tmp/{sample}_R2_Trimmingstats_after_trimming", 
        assem = "tmp/AssemblyQC.txt",
        mlst = "tmp/MLST.tsv"
    output:
        "stats/Summary_{sample}.tsv"
    params:
        sample = "{sample}"
    shell:
        "R1.before=$(cat {input.R1before}); R1.after=$(cat {input.R1after})"
        "&& R2.before=$(cat {input.R2before}); R2.after=$(cat {input.R2after})"
        "&& Total=$((R1.after+R2.after))"
        "&& assemstats=$(grep {params.sample} {input.assem}"
        "&& EstGenomeSize=$(cut -f 16,16 $assemstats)"
        "&& AvCoverage=$(($Total/$EstGenomeSize))"
        "&& GC=$(cut -f 17,17 $assemstats);N50=$(cut -f 18,18 $assemstats)"
        "&& MLST=$(grep {params.sample} {input.mlst}| cut -f 1,1 --complement)"
        "for i in $R1.before $R1.after $R2.before $R2.after $Total $EstGenomeSize $AvCoverage $GC $N50 $MLST; do echo -en $i '\t'>> {output}; done"
        "echo  >> {output}"

rule sumstat:
     input:
        expand("stats/Summary_{sample}.tsv", sample=SAMPLES)
     output:
        "stats/Summary.tsv"
     shell:
        "echo -e 'Reads\t#R1-Bases before trimming \t#R1-Bases after trimming\t#R2-Bases before trimming \t#R2-Bases after trimming\t#Bases Total\tEst. Genome size\tAv. coverage\t%GC\tN50\tspecies (MLST)\tMLST' > {output} "
        "&& cat {input} >> {output} "