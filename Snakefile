
SAMPLES = ["ENH-JSC-EDH-100490"]
R = ["R1","R2"]
SEQTKPARAM = ["-q 0.01"]


rule all:
    input:
        "stats/Trimmingstats.tsv"

rule fastqc:
    input:
        "{prefix}.fastq.gz"
    output:
        "{prefix}_fastqc.zip"
    shell:
        "fastqc {input}"


rule trim:
     input:
        "data/{sample}_{r}.fastq.gz"
     output:
        temp("trimmed/{sample}_{r}.fastq.gz")
     params:
        "-q 0.01"
     shell: 
        "seqtk trimfq {params} {input} > {output}"

rule trimstat:
     input:
        files=expand('trimmed/{sample}_{r}.fastq.gz', sample=SAMPLES, r=R)
     output:
       "stats/Trimmingstats.tsv"
     run:
        with open(output[0], 'w') as out:
            for i in input:
                out.write(i)
                out.write("\t")
#               sample = i.split('.')[0]
#                for line in open(i):
#                    out.write(sample + ' ' + line)

