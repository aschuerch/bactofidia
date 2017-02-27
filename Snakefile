import os
import random
from snakemake.utils import min_version

min_version("3.2")
configfile: "config.yaml"

#Minimum length of contigs
minlen = config["minlen"]
#Minimum coverage depth of contigs
mincov = config["mincov"]

seqtkversion = config["seqtk"]["version"] 
seqtkparams = config["seqtk"]["params"]

#Assembly with SPAdes
spadesparams = config["SPAdes"]["params"]
spadesversion = config["SPAdes"]["version"]
krange = config.get("krange")
#krange = "57,97,127" #adjust if you would like to select a different kmer range for spades assembly
versiontag = config["SPAdes"]["versiontag"]

# Prokka
prokkaparams = config["prokka"]["params"] 
prokkaversion = config["prokka"]["version"]
# mlst
mlstparams = config["mlst"]["params"]
mlstversion = config["mlst"]["version"]

# QUAST  
quastparams = config["QUAST"]["params"]
quastversion = config["QUAST"]["version"]

# Kraken
krakenparams = config["kraken"]["params"] 
krakenversion = config["kraken"]["version"] 

# Krona
kronaversion = config["krona"]["version"]
# Checkm
checkmparams = config["checkm"]["params"]
checkmversion = config["checkm"]["version"]

rand=random.random()

# python 2 virtual environments
os.system ( "conda create -y -n bactofidia{}  python=2.7  checkm-genome={} quast={} seqtk={} spades={} prokka={} mlst={} kraken={} krona={}".format ( rand, checkmversion, quastversion, seqtkversion, spadesversion, prokkaversion, mlstversion, krakenversion, kronaversion)    )

os.system ( "conda-env export -n bactofidia{} > condaenv{}.yml".format (rand,rand)) 


# Collect samples

SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")
SAMPLES = set(SAMPLES)


onsuccess:
    # delete virtual environment
    os.system ( "conda-env remove -n bactofidia{}".format (rand)) 
    print("Workflow finished!")


onerror:
    # delete virtual environment
    os.system ( "conda-env remove -n bactofidia{}".format (rand)) 
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
    params:
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"
        "source deactivate"

rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("trimmed/{sample}_{r}.fastq")
     params:
        p=seqtkparams,
        rand=rand
     shell: 
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "seqtk trimfq {params.p} {input} > {output}"
        "source deactivate"

rule fastqc_after:
    input:
        "trimmed/{sample}_{r}.fastq"
    output:
         temp("stats/{sample}_{r}_Trimmingstats_after_trimming")
    params:
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"
        "source deactivate"

rule trimstat:
    input:
        after=expand("stats/{sample}_{r}_Trimmingstats_after_trimming", sample=SAMPLES, r=R),
        before=expand("stats/{sample}_{r}_Trimmingstats_before_trimming", sample=SAMPLES, r=R)
    output:
        "stats/Trimmingstats.tsv"
    params:
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "echo -e 'Reads\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high' > {output} ;"
        "cat {input.before} >> {output} ;"
        "cat {input.after} >> {output} ;"
        "source deactivate"
       
rule spades:
    input: 
        R1="trimmed/{sample}_R1.fastq",
        R2="trimmed/{sample}_R2.fastq"
    output:
        "assembly/{sample}/scaffolds.fasta"
    params:
        rand=rand,
        spadesparams = spadesparams,
        kmer = krange,
        cov = mincov,
        spadesversion = spadesversion,
        outfolder = "assembly/{sample}"
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "python {params.spadesversion} -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"
        "source deactivate"

rule rename:
    input:
        "assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna"
    params:
        rand=rand,
        minlen = minlen,
        versiontag = "{sample}:"+versiontag
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.versiontag}/g > {output}"
        "source deactivate"

rule annotation:
   input:
        "scaffolds/{sample}.fna"
   output:
        "annotation/{sample}.gff"
   params:
        rand=rand,
        dir = "annotation",
        params = prokkaparams,
        prefix = "{sample}",
   shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "prokka --force --prefix {params.prefix} --outdir {params.dir} {params.params} {input} "
        "source deactivate"
        
rule taxonomy_1:
   input:
        "scaffolds/{sample}.fna"
   output:
        temp("taxonomy/{sample}.krakenout")
   params:
        p=krakenparams,
        rand=rand
   shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "kraken {params.p} --output {output} --fasta_input {input}"
        "source deactivate"

rule taxonomy_2:
    input:
        "taxonomy/{sample}.krakenout"
    output:
        temp("taxonomy/{sample}.kronain")
    params: 
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "cut -f2,3 {input} > {output}"
        "source deactivate"

rule taxonomy_3:
    input:
        "taxonomy/{sample}.kronain"
    output:
        "stats/Taxonomy_{sample}.html"
    params: 
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "ktImportTaxonomy {input} -o {output}"
        "source deactivate"

rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    params:
        mlstparams,
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "mlst {params} {input} >> {output}"
        "source deactivate"

rule quast:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        html = "stats/AssemblyQC.html",
        tsv = "stats/AssemblyQC.tsv"
    params:
        rand=rand,
        outfolder = "stats/quasttemp",
        minlen = minlen,
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "quast.py --min-contig {params.minlen} -o {params.outfolder} {input}"
        " && mv stats/quasttemp/report.html {output.html}"
        " && mv stats/quasttemp/report.tsv {output.tsv}"
        " && rm -r stats/quasttemp"
        "source deactivate"
 
rule resfinder:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES) 
    output:
        "stats/ResFinder.tsv"
    params: 
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "abricate {input} > {output}" 
        "&& sed -i 's/scaffolds\///g' {output}"
        "source deactivate"

rule checkm:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        file = "stats/SpeciesDetermination.tsv",
        folder = temp("checkm")
    params:
        checkmparams,
        rand=rand
    shell:
        "set +u; source activate bactofidia{params.rand}; set -u;"
        "checkm lineage_wf scaffolds {output.folder} {params} --tab_table -x fna > {output.file} ;"
        "source deactivate"

