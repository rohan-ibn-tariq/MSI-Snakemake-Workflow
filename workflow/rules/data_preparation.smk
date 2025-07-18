"""
# Data Preparation & Related Rules
# =========================
# Handles genome acquisition, repeat discovery, and MSI region identification
"""


rule get_genome:
    output:
        "workflow/resources/reference-sequence/chr{chrom}/genome.fasta"
    conda:
        "../envs/get_genome.yaml"
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        chromosome=lambda wildcards: [] if wildcards.chrom == "all" else [wildcards.chrom]
    log:
        "workflow/logs/data-preparation/get_genome_chr{chrom}.log"
    wrapper:
        "v6.2.0/bio/reference/ensembl-sequence"


rule bwa_index:
    input:
        "workflow/resources/reference-sequence/chr{chrom}/genome.fasta"
    output:
        "workflow/resources/reference-sequence/chr{chrom}/genome.fasta.bwt"
    conda:
        "../envs/simulation_and_variation.yaml"
    log:
        "workflow/logs/data-preparation/bwa_index_chr{chrom}.log"
    shell:
        "bwa index {input} > {log} 2>&1"


rule run_pytrf:
    input:
        fasta="workflow/resources/reference-sequence/chr{chrom}/genome.fasta"
    output:
        csv="results/data-preparation/repeats/chr{chrom}_pytrf_output.csv"
    log:
        "workflow/logs/data-preparation/run_pytrf_chr{chrom}.log"
    conda:
        "../envs/pytrf.yaml"
    params:
        min_repeats=config["msi_analysis"]["min_repeats"],
        fmt="csv"
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
        csv="results/data-preparation/repeats/chr{chrom}_pytrf_output.csv"
    output:
        bed="results/data-preparation/ms-bed/chr{chrom}_sample_test.bed"
    log:
        "workflow/logs/data-preparation/csv_to_bed_chr{chrom}.log"
    params:
        script="workflow/scripts/csv2bed_ucsc.py"
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        "python {params.script} --csv {input.csv} --bed {output.bed} > {log} 2>&1"


rule get_complete_annotation:
    output:
        "workflow/resources/annotation/complete.gtf"
    params:
        species=config["reference"]["species"],
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        fmt="gtf",
    log:
        "workflow/logs/data-preparation/get_complete_annotation.log"
    wrapper:
        "v6.2.0/bio/reference/ensembl-annotation"


rule split_annotation_by_chromosome:
    input:
        "workflow/resources/annotation/complete.gtf"
    output:
        "workflow/resources/annotation/chr{chrom}.gtf"
    log:
        "workflow/logs/data-preparation/split_annotation_chr{chrom}.log"
    shell:
        """
        if [ "{wildcards.chrom}" = "all" ]; then
            cp {input} {output} 2> {log}
        else
            grep "^{wildcards.chrom}[[:space:]]" {input} > {output} 2> {log}
        fi
        """
