SAMPLES,R, = glob_wildcards("data/{id}_{r}.fastq.gz")

#SAMPLES = ["ENH-JSC-EDH-100490"]
#R = ["R1","R2"]

seqtkparam = "-q 0.01  -l 30 -b 0 -e 0"
spadesparam = "--only-assembler --careful --threads 8"
spadesversion = "3.6.2"
minlen = 500
mincov = 10
krange = "57,97,127"
quastparam = "--min-contig 500"
version = "V022017SPAdes"+spadesversion

prokkaparam = "--centre UMCU --compliant"
krakenparam = "--quick --db /hpc/local/CentOS7/dla_mm/tools/kraken/db/minikraken_20140330/"
checkmparam = ""

def determine_spadespath(spadesversion):
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
       expand ("taxonomy/{sample}.html", sample=SAMPLES)


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
        seqtkparam
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
        "grep "" {input.before} >> {output} ;"
        "grep "" {input.after} >> {output} ;"
       
rule spades:
    input: 
        R1="trimmed/{sample}_R1.fastq",
        R2="trimmed/{sample}_R2.fastq"
    output:
        "assembly/{sample}/scaffolds.fasta"
    params:
        spadesparam = spadesparam,
        kmer = krange,
        cov = mincov,
        spadesversion = determine_spadespath(spadesversion),
        outfolder = "assembly/{sample}"
    shell:
        "python {params.spadesversion} -1 {input.R1} -2 {input.R2} -o {params.outfolder} -k {params.kmer} --cov-cutoff {params.cov} {params.spadesparam}"

rule rename:
    input:
        "assembly/{sample}/scaffolds.fasta"
    output:
        "scaffolds/{sample}.fna"
    params:
        minlen = minlen,
        version = "{sample}:"+version
    shell:
        "seqtk seq -L {params.minlen} {input} | sed  s/NODE/{params.version}/g > {output}"

rule annotation:
   input:
        "scaffolds/{sample}.fna"
   output:
        "annotation/{sample}.gff"
   params:
        dir = "annotation",
        param = prokkaparam,
        prefix = "{sample}"
   shell:
        "prokka --force --prefix {params.prefix} --outdir {params.dir} {params.param} {input} "
        
rule taxonomy_1:
   input:
        "scaffolds/{sample}.fna"
   output:
        temp("taxonomy/{sample}.krakenout")
   params:
        krakenparam
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
        "taxonomy/{sample}.html"
   shell:
        "ktImportTaxonomy {input} -o {output}"


rule mlst:
    input:
        expand("scaffolds/{sample}.fna", sample=SAMPLES)
    output:
        "stats/MLST.tsv"
    shell:
        "mlst --quiet --nopath {input} >> {output}"

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
        "&& sed -i 's/scaffolds//g' {output}"


rule checkfinder:
    input:
       "scaffolds/"
    output:
       file = "stats/SpeciesDetermination.tsv",
       folder = temp("checkm")
    shell:
       "set +u; source activate checkm; set -u;"
       "checkm lineage_wf {input} {output.folder} -x fna > {output.file} ;"
       "source deactivate"

