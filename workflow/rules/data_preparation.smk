"""
# Data Preparation & Related Rules
# =========================
# Handles genome acquisition, repeat discovery, and MSI region identification
"""


rule get_genome:
    output:
        "workflow/resources/reference-sequence/genome.fasta",
    conda:
        "../envs/get_genome.yaml"
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="112",
        chromosome=["22"],
    log:
        "workflow/logs/data-preparation/get_genome.log",
    wrapper:
        "v6.2.0/bio/reference/ensembl-sequence"


rule bwa_index:
    input:
        "workflow/resources/reference-sequence/genome.fasta",
    output:
        "workflow/resources/reference-sequence/genome.fasta.bwt",
    conda:
        "../envs/simulation_and_variation.yaml"
    log:
        "workflow/logs/data-preparation/bwa_index.log",
    shell:
        "bwa index {input} > {log} 2>&1"


rule run_pytrf:
    input:
        fasta="workflow/resources/reference-sequence/genome.fasta",
    output:
        csv="results/data-preparation/repeats/pytrf_output.csv",
    log:
        "workflow/logs/data-preparation/run_pytrf.log",
    conda:
        "../envs/pytrf.yaml"
    params:
        min_repeats=(5, 5, 5, 5, 5, 5),  # mono, di, tri, tetra, penta, hexa
        fmt="csv",
    shell:
        """
        pytrf findstr \
            -r {params.min_repeats[0]} {params.min_repeats[1]} {params.min_repeats[2]} \
               {params.min_repeats[3]} {params.min_repeats[4]} {params.min_repeats[5]} \
            -f {params.fmt} \
            -o {output.csv} \
            {input.fasta} > {log} 2>&1
        """


rule csv_to_bed:
    input:
        csv="results/data-preparation/repeats/pytrf_output.csv",
    output:
        bed="results/data-preparation/ms-bed/genome.bed",
    log:
        "workflow/logs/data-preparation/csv_to_bed.log",
    params:
        script="workflow/scripts/csv2bed_ucsc.py",
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        "python {params.script} --csv {input.csv} --bed {output.bed} > {log} 2>&1"


rule create_test_bed:
    input:
        "results/data-preparation/ms-bed/genome.bed",
    output:
        "results/data-preparation/ms-bed/sample_test.bed",
    log:
        "workflow/logs/data-preparation/create_test_bed.log",
    params:
        num_entries=4000,
    shell:
        "head -n {params.num_entries} {input} > {output} 2> {log}"
