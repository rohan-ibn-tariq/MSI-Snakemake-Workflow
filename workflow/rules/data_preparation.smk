"""
# Data Preparation & Related Rules
# =========================
# Handles genome acquisition, repeat discovery, and MSI region identification
"""


rule get_genome:
    output:
        "workflow/resources/genome.fasta",
    conda:
        "../envs/get_genome.yaml"
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="112",
        chromosome=["22"],
    log:
        "workflow/logs/get_genome.log",
    wrapper:
        "v6.2.0/bio/reference/ensembl-sequence"


rule run_pytrf:
    input:
        fasta="workflow/resources/genome.fasta"
    output:
        csv="results/repeats/pytrf_output.csv"
    log:
        "workflow/logs/run_pytrf.log"
    conda:
        "../envs/pytrf.yaml"
    params:
        min_repeats=(5, 5, 5, 5, 5, 5),  # mono, di, tri, tetra, penta, hexa
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
        csv="results/repeats/pytrf_output.csv"
    output:
        bed="results/ms-bed/genome.bed"
    log:
        "workflow/logs/csv_to_bed.log"
    params:
        script="workflow/scripts/csv2bed_ucsc.py"
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        "python {params.script} --csv {input.csv} --bed {output.bed} > {log} 2>&1"


rule create_test_bed:
    input: "results/ms-bed/genome.bed"
    output: "results/ms-bed/sample_test.bed"
    log:
        "workflow/logs/create_test_bed.log"
    params: num_entries=1500
    shell: "head -n {params.num_entries} {input} > {output} 2> {log}"
