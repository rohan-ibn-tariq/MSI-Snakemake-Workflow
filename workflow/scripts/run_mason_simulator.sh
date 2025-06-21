#!/bin/bash
set -euo pipefail

# Set file descriptor limit
# NOTE: This is pysam workaround, 
# see: https://github.com/seqan/seqan/issues/2490#issuecomment-1893933011
ulimit -n "${snakemake_params[file_limit]}"

# Run mason simulator
mason_simulator \
    -ir "${snakemake_input[fasta]}" \
    -iv "${snakemake_input[vcf]}" \
    -n "${snakemake_params[num_fragments]}" \
    --illumina-read-length "${snakemake_params[read_length]}" \
    --fragment-mean-size "${snakemake_params[fragment_mean_size]}" \
    --fragment-size-std-dev "${snakemake_params[fragment_size_stddev]}" \
    --seed "${snakemake_params[seed]}" \
    --num-threads "${snakemake[threads]}" \
    -o "${snakemake_output[fastq_left]}" \
    -or "${snakemake_output[fastq_right]}" \
    -oa "${snakemake_output[bam]}"

# Sort and index BAM file
samtools sort -@ "${snakemake[threads]}" -o "${snakemake_output[bam]}.tmp" "${snakemake_output[bam]}"
mv "${snakemake_output[bam]}.tmp" "${snakemake_output[bam]}"

samtools index "${snakemake_output[bam]}"
