"""
# Simulation and Processing of Variants for MSI Detection Related Rules
# =======================================================================
# This Snakemake workflow simulates variants using Mason, processes them to inject
# microsatellite instability (MSI) variants, and prepares them for further analysis.
# It includes rules for variant simulation, VCF format correction, and read simulation.
# The workflow is designed to be modular and can be extended with additional processing steps.
"""


rule simulate_variants_mason:
    input:
        ref="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
    output:
        vcf="results/simulation-processing/simulated-mason-variant/chr{chrom}_sim_raw.vcf",
    params:
        seed=50,
        haplotype=2,
        small_indel_rate=0.001,
    log:
        "workflow/logs/simulation-processing/simulate_variants_mason_chr{chrom}.log",
    conda:
        "../envs/simulation_and_variation.yaml"
    shell:
        """
        mason_variator \
            -ir {input.ref} \
            -ov {output.vcf} \
            --seed {params.seed} \
            -n {params.haplotype} \
            --small-indel-rate {params.small_indel_rate} \
            > {log} 2>&1
        """


rule fix_mason_vcf_format:
    input:
        raw_vcf="results/simulation-processing/simulated-mason-variant/chr{chrom}_sim_raw.vcf",
    output:
        vcf="results/simulation-processing/simulated-mason-variant/chr{chrom}_sim.vcf",
    params:
        format_tag="GT",
        description="Genotype",
    log:
        "workflow/logs/simulation-processing/fix_mason_vcf_format_chr{chrom}.log",
    shell:
        """
        awk -v format_tag="{params.format_tag}" \
            -v description="{params.description}" \
            -f workflow/scripts/fix_mason_vcf.awk \
            {input.raw_vcf} > {output.vcf} 2> {log}
        """


rule process_indels:
    input:
        bed="results/data-preparation/ms-bed/chr{chrom}_sample_test.bed",
        vcf="results/simulation-processing/simulated-mason-variant/chr{chrom}_sim.vcf",
        ref_fasta="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
    output:
        "results/simulation-processing/msi-indels-vcf/chr{chrom}_sim_indels.vcf",
    log:
        "workflow/logs/simulation-processing/process_indels_chr{chrom}.log",
    conda:
        "../envs/data_processing_and_injection.yaml"
    params:
        ins_rate=0.30,
        del_rate=0.30,
        boost_rate=0.70,
        seed=50,
    shell:
        """
        PYTHONHASHSEED=0 python workflow/scripts/extending_vcf_for_indels.py \
            --bed {input.bed} \
            --vcf-in {input.vcf} \
            --vcf-out {output} \
            --ref-fasta {input.ref_fasta} \
            --ins-rate {params.ins_rate} \
            --del-rate {params.del_rate} \
            --boost-rate {params.boost_rate} \
            --seed {params.seed} \
            --verbose > {log} 2>&1
        """


# NOTE: KEEP THIS COMMENT BLOCK
# Filter out all BND's if sorting active in process indels
# as they are not required for MSI detection.
# NOTE: In case required for more realistic simulation,
# the sorting in process indels
# should not break the BND chunks in 6.
# rule filter_vcf_for_mason:
#     input:
#         "results/msi-indels-vcf/genome_sim_indels.vcf"
#     output:
#         "results/msi-indels-vcf/genome_filtered.vcf"
#     log:
#         "workflow/logs/filter_vcf_for_mason.log"
#     conda:
#         "../envs/data_processing_and_injection.yaml"
#     shell:
#         """
#         echo "Before filtering: $(grep -v '^#' {input} | wc -l) variants" > {log}
#         bcftools view -e 'INFO/SVTYPE="BND"' -Ov -o {output} {input} 2>> {log}
#         echo "After filtering: $(grep -v '^#' {output} | wc -l) variants" >> {log}
#         """


rule simulate_reads_mason:
    input:
        fasta="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
        vcf="results/simulation-processing/msi-indels-vcf/chr{chrom}_sim_indels.vcf",
    output:
        fastq_left="results/simulation-processing/simulated-fq/chr{chrom}_left_reads.fq",
        fastq_right="results/simulation-processing/simulated-fq/chr{chrom}_right_reads.fq",
    params:
        num_fragments=config["simulation"]["num_fragments"],
        read_length=config["simulation"]["read_length"],
        fragment_mean_size=350,
        fragment_size_stddev=50,
        seed=config["simulation"]["seed"],
        file_limit=65536,
    threads: 4
    log:
        "workflow/logs/simulation-processing/simulate_reads_mason_chr{chrom}.log",
    conda:
        "../envs/simulation_and_variation.yaml"
    shell:
        """
        # Set file descriptor limit
        # NOTE: This is pysam workaround, 
        # see: https://github.com/seqan/seqan/issues/2490#issuecomment-1893933011
        ulimit -n {params.file_limit}
        
        mason_simulator \
            -ir {input.fasta} \
            -iv {input.vcf} \
            -n {params.num_fragments} \
            --illumina-read-length {params.read_length} \
            --fragment-mean-size {params.fragment_mean_size} \
            --fragment-size-std-dev {params.fragment_size_stddev} \
            --seed {params.seed} \
            --num-threads {threads} \
            -o {output.fastq_left} \
            -or {output.fastq_right} > {log} 2>&1
        """


rule align_with_bwa:
    input:
        fastq_left="results/simulation-processing/simulated-fq/chr{chrom}_left_reads.fq",
        fastq_right="results/simulation-processing/simulated-fq/chr{chrom}_right_reads.fq",
        ref="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
        index="workflow/resources/reference-sequence/chr{chrom}/genome.fasta.bwt",
    output:
        bam="results/simulation-processing/simulated-bam/chr{chrom}_alignments.bam",
    threads: 4
    log:
        "workflow/logs/simulation-processing/align_with_bwa_chr{chrom}.log",
    conda:
        "../envs/simulation_and_variation.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fastq_left} {input.fastq_right} | \
        samtools sort -@ {threads} -o {output.bam} > {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """


rule use_msi_variants_as_candidates:
    input:
        "results/simulation-processing/msi-indels-vcf/chr{chrom}_sim_indels.vcf",
    output:
        "results/simulation-processing/candidates/chr{chrom}_candidates.vcf",
    log:
        "workflow/logs/simulation-processing/use_msi_variants_as_candidates_chr{chrom}.log",
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        """
        bcftools sort {input} -Ov -o {output} > {log} 2>&1
        """


rule estimate_alignment_properties:
    input:
        bam="results/simulation-processing/simulated-bam/chr{chrom}_alignments.bam",
        ref="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
    output:
        "results/simulation-processing/alignment-properties/chr{chrom}_sample.alignment-properties.json",
    conda:
        "../envs/varlociraptor.yaml"
    log:
        "workflow/logs/simulation-processing/estimate_alignment_properties_chr{chrom}.log",
    shell:
        """
        varlociraptor estimate alignment-properties \
            {input.ref} \
            --bams {input.bam} > {output} 2> {log}
        """


rule varlociraptor_preprocess:
    input:
        bam="results/simulation-processing/simulated-bam/chr{chrom}_alignments.bam",
        ref="workflow/resources/reference-sequence/chr{chrom}/genome.fasta",
        candidates="results/simulation-processing/candidates/chr{chrom}_candidates.vcf",
        alignment_properties="results/simulation-processing/alignment-properties/chr{chrom}_sample.alignment-properties.json",
    output:
        "results/simulation-processing/varlociraptor-preprocess/chr{chrom}_observations.bcf",
    threads: 4
    conda:
        "../envs/varlociraptor.yaml"
    log:
        "workflow/logs/simulation-processing/varlociraptor_preprocess_chr{chrom}.log",
    shell:
        """
        varlociraptor preprocess variants {input.ref} \
            --candidates {input.candidates} \
            --bam {input.bam} \
            --alignment-properties {input.alignment_properties} \
            --atomic-candidate-variants \
            --pairhmm-mode homopolymer \
            > {output} \
            2> {log}
        """


rule call_variants_varlociraptor:
    input:
        obs="results/simulation-processing/varlociraptor-preprocess/chr{chrom}_observations.bcf",
        scenario="workflow/resources/scenario.yaml",
    output:
        "results/simulation-processing/varlociraptor-calls/chr{chrom}_variants.bcf",
    threads: 4
    log:
        "workflow/logs/simulation-processing/call_variants_varlociraptor_chr{chrom}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        """
        varlociraptor call variants \
            --output {output} \
            --omit-homopolymer-artifact-detection \
            --omit-alt-locus-bias \
            generic \
            --scenario {input.scenario} \
            --obs sample={input.obs} \
            2> {log}
        """
