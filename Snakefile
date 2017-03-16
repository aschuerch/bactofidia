import os
import sys
import glob
import subprocess
import shutil

from snakemake.utils import min_version

min_version("3.9")
configfile: "config.yaml"

##to build a virt env
seqtkversion = config["seqtk"]["version"] 
spadesversion = config["SPAdes"]["version"]
versiontag = "V022017SPAdes"+spadesversion


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

# abricate
abricateversion = config["abricate"]["version"]

def kmer_determination():
    if (config.get("krange")):
        kmer = config.get("krange")
    else:
        kmer = config["SPAdes"] ["krange"]
    return kmer

# virtual environments
virtenvs = "checkm-genome={} quast={} seqtk={} spades={} prokka={} mlst={} kraken={} krona={} abricate={}".format ( checkmversion, quastversion, seqtkversion, spadesversion, prokkaversion, mlstversion, krakenversion, kronaversion, abricateversion).split()

os.makedirs ("virtenvs", exist_ok=True )

def spec_virtenv(program):
    """ Create conda enviroment for every job. """
    if shutil.which("conda") is None:
        raise CreateCondaEnvironmentException("The 'conda' command is not available in $PATH.")

    i = (list( filter(lambda x: program in x, virtenvs))[0])         
   # stdout = open("virtenvs/{}.txt".format(program),"wb")

    print (program, "version", i)
    try: 
  #      x = subprocess.check_output(["conda", "env", "export", "--prefix", "virtenvs/"+program]) 
   #     stdout.write(x)
   # except: 
        subprocess.check_output([ "conda", "create","-y", "--prefix", "virtenvs/"+program, i, "--copy"]  , stderr=subprocess.STDOUT)
    #    x = subprocess.check_output(["conda", "env", "export", "--prefix", "virtenvs/"+program]) 
     #   stdout.write(x)
    except: 
        print ('Virtual environment already present')

    return "virtenvs/{}".format(program)




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
       "stats/Trimmingstats.tsv",
       "stats/MLST.tsv",
       "stats/AssemblyQC.html",
       "stats/AssemblyQC.tsv",
       "stats/ResFinder.tsv",
       "stats/SpeciesDetermination.tsv",
       expand ("scaffolds/{sample}.fna", sample=SAMPLES),
       expand ("annotation/{sample}.gff", sample=SAMPLES),
#       expand ("stats/Taxonomy_{sample}.html", sample=SAMPLES)


rule fastqc_before:
    input:
        "data/{sample}_{r}.fastq.gz"
    output:
        temp("stats/{sample}_{r}_Trimmingstats_before_trimming")
    params:
        virtenv=spec_virtenv("seqtk")
    shell:
        "set +u; source activate {params.virtenv}; set -u;"
        "&& seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"
        "&& set +u; source deactivate; set -u;"

rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("trimmed/{sample}_{r}.fastq")
     params:
        p = config["seqtk"]["params"],
        virtenv=spec_virtenv("seqtk")
     shell: 
        "set +u; source activate {params.virtenv}; set -u;"
        "&& seqtk trimfq {params.p} {input} > {output} "
        "&& set +u; source deactivate; set -u;"

rule fastqc_after:
    input:
        "trimmed/{sample}_{r}.fastq"
    output:
         temp("stats/{sample}_{r}_Trimmingstats_after_trimming")
    params:
        virtenv=spec_virtenv("seqtk")
    shell:
        "set +u; source activate {params.virtenv}; set -u;"
        "&& seqtk fqchk {input} | grep ALL | sed 's/ALL//g' >> {output}"
        "&& set +u; source deactivate; set -u;"       

rule trimstat:
    input:
        after=expand("stats/{sample}_{r}_Trimmingstats_after_trimming", sample=SAMPLES, r=R),
        before=expand("stats/{sample}_{r}_Trimmingstats_before_trimming", sample=SAMPLES, r=R)
    output:
        "stats/Trimmingstats.tsv"
    shell:
        "echo -e 'Reads\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high' > {output} "
        "&& echo -en {input.before}  >> {output} "
        "&& cat {input.before} >> {output} "
        "&& echo -en {input.after} >> {output} "
        "&& cat {input.after} >> {output} "
       
rule spades:
    input: 
        R1="trimmed/{sample}_R1.fastq",
        R2="trimmed/{sample}_R2.fastq"
    output:
        "assembly/{sample}/scaffolds.fasta"
    params:
        spadesparams = config["SPAdes"]["params"],
        kmer = kmer_determination(),
        cov = config["mincov"],
        outfolder = "assembly/{sample}",
        virtenv=spec_virtenv("spades")
    shell:
        "set +u; source activate {params.virtenv}; set -u;"
        "&& spades.py -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparams}"
        "&& set +u; source deactivate; set -u;"

rule rename:
    input:
        "assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna",
    params:
        minlen = config["minlen"],
        versiontag = "{sample}_"+versiontag,
        virtenv=spec_virtenv("seqtk")
    shell:
        "set +u; source activate {params.virtenv}; set -u;"
        "&& seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.versiontag}/g > {output}"
        "&& set +u; source deactivate; set -u;"
      
rule annotation:
    input:
        "scaffolds/{sample}.fna"
    output:
        "annotation/{sample}.gff",
    params:
        dir = "annotation",
        params = config["prokka"]["params"],
        prefix = "{sample}",
        virtenv=spec_virtenv("prokka")
    shell:
        "set +u; source activate {params.virtenv}; set -u;"
        "&& prokka --force --prefix {params.prefix} --outdir {params.dir} {params.params} {input} "
        "&& set +u; source deactivate; set -u;"
        
rule taxonomy_1:
    input:
        "scaffolds/{sample}.fna"
    output:
        temp("taxonomy/{sample}.krakenout")
    params:
        p=config["kraken"]["params"],
        virtenv=spec_virtenv("kraken")
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& kraken {params.p} --output {output} --fasta_input {input}"
        "&& set +u; source deactivate; set -u;"

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
    params:
        virtenv=spec_virtenv("krona")
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& ktImportTaxonomy {input} -o {output}"
        "&& set +u; source deactivate; set -u;"


rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    params:
        config["mlst"]["params"],
        virtenv=spec_virtenv("mlst")
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& mlst {params} {input} >> {output}"
        "&& set +u; source deactivate; set -u;"

rule quast:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        html = "stats/AssemblyQC.html",
        tsv = "stats/AssemblyQC.tsv"
    params:
        outfolder = "stats/quasttemp",
        p = config["QUAST"]["params"],
        virtenv=spec_virtenv("quast")
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& quast {params.p} -o {params.outfolder} {input}"
        "&& mv stats/quasttemp/report.html {output.html}"
        "&& mv stats/quasttemp/report.tsv {output.tsv}"
        "&& rm -r stats/quasttemp"
        "&& set +u; source deactivate; set -u;"

rule resfinder:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES) 
    output:
        "stats/ResFinder.tsv"
    params:
        virtenv = spec_virtenv("abricate")
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& abricate {input} > {output}" 
        "&& sed -i 's/scaffolds\///g' {output}"
        "&& set +u; source deactivate; set -u;"

rule checkm:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        file = "stats/SpeciesDetermination.tsv",
        folder = temp("checkm")
    params:
        config["checkm"]["params"],
        virtenv = spec_virtenv("checkm")        
    shell:
        "set +u; source activate  {params.virtenv}; set -u;"
        "&& checkm lineage_wf scaffolds {output.folder} {params} --tab_table -x fna > {output.file}"
        "&& set +u; source deactivate; set -u;"
