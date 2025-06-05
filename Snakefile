rule all:
    input:
        "resources/genome.fasta"


rule get_genome:
    output:
        "resources/genome.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="112",
        chromosome=["22"],
    log:
        "logs/get_genome.log",
    wrapper:
        "v6.2.0/bio/reference/ensembl-sequence"


rule dat_to_bed:
    input:
        dat="results/trf/{sample}.dat"
    output:
        bed="results/msi/{sample}.bed"
    log:
        "logs/{sample}.log"
    params:
        script="scripts/dat2bed_ucsc.py"
    conda:
        "environment/bed.yaml"
    shell:
        "python {params.script} --dat {input.dat} --bed {output.bed} > {log} 2>&1"
