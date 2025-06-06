rule all:
    input:
        [
            "resources/genome.fasta",
            "results/simulated/genome_sim.vcf",
            "genome.fasta.2.7.7.80.10.50.500.dat",
            "results/msi/genome.bed"
        ]


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


rule simulate_variants_mason:
    input:
        ref="resources/genome.fasta"
    output:
        vcf="results/simulated/genome_sim.vcf"
    params:
        seed=50,
        num_variants=1000
    log:
        "logs/simulate_variants.log"
    shell:
        """
        mason_variator \
            -ir {input.ref} \
            -ov {output.vcf} \
            --seed {params.seed} \
            -n {params.num_variants} \
            --small-indel-rate 0.001        """


rule run_trf_on_reference:
    input:
        ref="resources/genome.fasta"
    output:
        dat="genome.fasta.2.7.7.80.10.50.500.dat"
    params:
        match=2,
        mismatch=7,
        delta=7,
        pm=80,
        pi=10,
        minscore=50,
        maxperiod=500,
        extra="-h",
    log:
        "logs/trf/trf.log"
    shell:
        """
        trf \
            {input.ref} \
            {params.match} \
            {params.mismatch} \
            {params.delta} \
            {params.pm} \
            {params.pi} \
            {params.minscore} \
            {params.maxperiod} \
            {params.extra} \
            > {log} 2>&1
        """


rule dat_to_bed:
    input:
        dat="genome.fasta.2.7.7.80.10.50.500.dat"
    output:
        bed="results/msi/genome.bed"
    log:
        "logs/msi.log"
    params:
        script="scripts/dat2bed_ucsc.py"
    shell:
        "python {params.script} --dat {input.dat} --bed {output.bed} > {log} 2>&1"
