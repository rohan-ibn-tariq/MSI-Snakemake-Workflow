"""
Snakemake rules for data analysis in the MSI detection workflow.
Handles .bcf to .vcf conversion and quantification of MSI variants.
"""


rule convert_bcf_to_vcf:
    input:
        "results/simulation-processing/varlociraptor-calls/chr{chrom}_variants.bcf",
    output:
        "results/data-analysis/varlociraptor-calls/chr{chrom}_variants.vcf",
    log:
        "workflow/logs/data-analysis/convert_bcf_to_vcf_chr{chrom}.log",
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        "bcftools view {input} -Ov -o {output} 2> {log}"


# rule convert_bcf_to_vcf:
#     input:
#         "results/simulation-processing/varlociraptor-calls/variants.bcf",
#     output:
#         "results/data-analysis/varlociraptor-calls/variants.vcf",
#     log:
#         "workflow/logs/simulation-processing/convert_bcf_to_vcf.log",
#     conda:
#         "../envs/data_processing_and_injection.yaml"
#     shell:
#         "bcftools view {input} -Ov -o {output} 2> {log}"


rule gtf_to_coding_bed:
    input:
        gtf="workflow/resources/annotation/chr{chrom}.gtf"
    output:
        bed="results/data-analysis/coding-regions/chr{chrom}_coding.bed"
    log:
        "workflow/logs/data-analysis/gtf_to_coding_bed_chr{chrom}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Extract CDS and add chr prefix (input is already chromosome-specific)
        awk '$3=="CDS" {{print "chr"$1"\\t"($4-1)"\\t"$5}}' {input.gtf} | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - > {output.bed} 2> {log}
        """


rule gtf_to_coding_bed_for_vcf:
    input:
        gtf="workflow/resources/annotation/chr{chrom}.gtf"
    output:
        bed="results/data-analysis/coding-regions/chr{chrom}_coding_for_vcf.bed"
    log:
        "workflow/logs/data-analysis/gtf_to_coding_bed_for_vcf_chr{chrom}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Remove chr prefix for VCF-BED intersections
        awk '$3=="CDS" {{print $1"\\t"($4-1)"\\t"$5}}' {input.gtf} | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge -i - > {output.bed} 2> {log}
        """


rule intersect_vcf_with_coding:
    input:
        vcf="results/data-analysis/varlociraptor-calls/chr{chrom}_variants.vcf",
        bed="results/data-analysis/coding-regions/chr{chrom}_coding_for_vcf.bed"
    output:
        vcf="results/data-analysis/coding-intersected/chr{chrom}_coding_variants.vcf"
    log:
        "workflow/logs/data-analysis/intersect_vcf_with_coding_chr{chrom}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.bed} -header > {output.vcf} 2> {log}
        """


rule intersect_msi_bed_with_coding:
    input:
        msi_bed="results/data-preparation/ms-bed/chr{chrom}_sample_test.bed",
        coding_bed="results/data-analysis/coding-regions/chr{chrom}_coding.bed"
    output:
        intersected_bed="results/data-analysis/coding-msi-bed/chr{chrom}_coding_msi.bed"
    log:
        "workflow/logs/data-analysis/intersect_msi_bed_with_coding_chr{chrom}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # Convert 5-column BED to 3-column BED for bedtools, then restore
        cut -f2,3,4 {input.msi_bed} > {output.intersected_bed}.temp3col
        bedtools intersect -a {output.intersected_bed}.temp3col -b {input.coding_bed} > {output.intersected_bed}.temp_intersected
        
        # Get line numbers of intersected regions
        bedtools intersect -a {output.intersected_bed}.temp3col -b {input.coding_bed} -c | \
        awk '$4>0 {{print NR}}' > {output.intersected_bed}.temp_lines
        
        # Extract original lines using line numbers
        while read line; do sed -n "${{line}}p" {input.msi_bed}; done < {output.intersected_bed}.temp_lines > {output.intersected_bed}
        
        # Cleanup
        rm {output.intersected_bed}.temp3col {output.intersected_bed}.temp_intersected {output.intersected_bed}.temp_lines
        """


rule quantify_msi_variants:
    input:
        vcf_varlociraptor="results/data-analysis/varlociraptor-calls/chr{chrom}_variants.vcf",
        vcf_ground_truth="results/simulation-processing/msi-indels-vcf/chr{chrom}_sim_indels.vcf",
    output:
        html_report="results/data-analysis/msi-quantification/chr{chrom}_report.html",
    log:
        "workflow/logs/data-analysis/quantify_msi_variants_chr{chrom}.log",
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        """
        python workflow/scripts/ms_instability_quantification.py \
            --ground-truth {input.vcf_ground_truth} \
            --varlociraptor-output {input.vcf_varlociraptor} \
            --html-output {output.html_report} \
            > {log} 2>&1
        """


# rule quantify_msi_variants:
#     input:
#         vcf_varlociraptor="results/data-analysis/varlociraptor-calls/variants.vcf",
#         vcf_ground_truth="results/simulation-processing/msi-indels-vcf/genome_sim_indels.vcf",
#     output:
#         html_report="results/data-analysis/msi-quantification/chr22_report.html",
#     log:
#         "workflow/logs/data-analysis/quantify_msi_variants.log",
#     conda:
#         "../envs/data_processing_and_injection.yaml"
#     shell:
#         """
#         python workflow/scripts/ms_instability_quantification.py \
#             --ground-truth {input.vcf_ground_truth} \
#             --varlociraptor-output {input.vcf_varlociraptor} \
#             --html-output {output.html_report} \
#             > {log} 2>&1
#         """


rule msi_analysis_full:
    input:
        vcf="results/data-analysis/varlociraptor-calls/chr{chrom}_variants.vcf",
        bed="results/data-preparation/ms-bed/chr{chrom}_sample_test.bed"
    output:
        results="results/data-analysis/msi-analysis-full/chr{chrom}/msi_results.json",
        html_report="results/data-analysis/msi-analysis-full/chr{chrom}/msi_report.html",
        debug_log="results/data-analysis/msi-analysis-full/chr{chrom}/msi_debug.log",
        quantification="results/data-analysis/msi-analysis-full/chr{chrom}/msi_quantification.json"
    threads: 4
    log:
        "workflow/logs/data-analysis/chr{chrom}_msi_analysis_full.log"
    benchmark:
        "benchmarks/data-analysis/chr{chrom}_msi_analysis_full.txt"
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        """
        cd workflow/scripts
        export PYTHONPATH="${{PYTHONPATH:-}}:$$(pwd)"
        python msi_quantification_module/__main__.py \
            --bed ../../{input.bed} \
            --vcf ../../{input.vcf} \
            --output ../../{output.results} \
            --html-report ../../{output.html_report} \
            --debug-log ../../{output.debug_log} \
            --quantification-output ../../{output.quantification} \
            --threads {threads} \
            > ../../{log} 2>&1
        """

rule msi_analysis_coding:
    input:
        vcf="results/data-analysis/coding-intersected/chr{chrom}_coding_variants.vcf",
        bed="results/data-analysis/coding-msi-bed/chr{chrom}_coding_msi.bed"
    output:
        results="results/data-analysis/msi-analysis-coding/chr{chrom}/msi_results.json",
        html_report="results/data-analysis/msi-analysis-coding/chr{chrom}/msi_report.html",
        debug_log="results/data-analysis/msi-analysis-coding/chr{chrom}/msi_debug.log",
        quantification="results/data-analysis/msi-analysis-coding/chr{chrom}/msi_quantification.json"
    threads: 4
    log:
        "workflow/logs/data-analysis/chr{chrom}_msi_analysis_coding.log"
    benchmark:
        "benchmarks/data-analysis/chr{chrom}_msi_analysis_coding.txt"
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        """
        cd workflow/scripts
        export PYTHONPATH="${{PYTHONPATH:-}}:$$(pwd)"
        python msi_quantification_module/__main__.py \
            --bed ../../{input.bed} \
            --vcf ../../{input.vcf} \
            --output ../../{output.results} \
            --html-report ../../{output.html_report} \
            --debug-log ../../{output.debug_log} \
            --quantification-output ../../{output.quantification} \
            --threads {threads} \
            > ../../{log} 2>&1
        """


# rule msi_analysis:
#     input:
#         vcf="results/data-analysis/varlociraptor-calls/variants.vcf",
#         bed="results/data-preparation/ms-bed/{sample}.bed"
#     output:
#         results="results/data-analysis/msi-analysis/{sample}_msi_results.json",
#         html_report="results/data-analysis/msi-analysis/{sample}_msi_report.html",
#         debug_log="results/data-analysis/msi-analysis/{sample}_msi_debug.log",
#         quantification="results/data-analysis/msi-analysis/{sample}_msi_quantification.json"
#     threads: 4
#     log:
#         "logs/data-analysis/{sample}_msi_analysis.log"
#     benchmark:
#         "benchmarks/data-analysis/{sample}_msi_analysis.txt"
#     shell:
#         """
#         # Set Python path to find the MSI module
#         export PYTHONPATH="${{PYTHONPATH}}:$(pwd)/workflow/scripts"
        
#         # Run MSI analysis with all outputs
#         python -m msi_quantification_module \
#             --bed {input.bed} \
#             --vcf {input.vcf} \
#             --output {output.results} \
#             --html-report {output.html_report} \
#             --debug-log {output.debug_log} \
#             --quantification-output {output.quantification} \
#             --threads {threads} \
#             > {log} 2>&1
#         """
