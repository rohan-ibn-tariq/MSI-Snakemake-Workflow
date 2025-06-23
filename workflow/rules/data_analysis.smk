"""
Snakemake rules for data analysis in the MSI detection workflow.
Handles .bcf to .vcf conversion and quantification of MSI variants.
"""


rule convert_bcf_to_vcf:
    input:
        "results/simulation-processing/varlociraptor-calls/variants.bcf",
    output:
        "results/data-analysis/varlociraptor-calls/variants.vcf",
    log:
        "workflow/logs/simulation-processing/convert_bcf_to_vcf.log",
    conda:
        "../envs/data_processing_and_injection.yaml"
    shell:
        "bcftools view {input} -Ov -o {output} 2> {log}"


rule quantify_msi_variants:
    input:
        vcf_varlociraptor="results/data-analysis/varlociraptor-calls/variants.vcf",
        vcf_ground_truth="results/simulation-processing/msi-indels-vcf/genome_sim_indels.vcf",
    output:
        html_report="results/data-analysis/msi-quantification/chr22_report.html",
    log:
        "workflow/logs/data-analysis/quantify_msi_variants.log",
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
