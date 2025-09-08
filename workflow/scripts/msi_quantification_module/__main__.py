"""
MSI Quantification Module - Main Entry Point

Command-line interface for MSI quantification analysis.
Usage: python -m msi_quantification_module [options]
"""


import argparse
import json
import os
from datetime import datetime


from msi_quantification_module.core import (
    load_vcf_variants,
    load_bed_regions, 
    bin_based_threaded_intersection,
    create_msi_quantification
)
from msi_quantification_module.debug import write_complete_debug_log, initialize_debug_log
from msi_quantification_module.dp_analysis import ( 
    prepare_variants_for_dp,
    run_af_evolution_analysis,
)
from msi_quantification_module.reports import (
    generate_msi_html_report, 
    generate_msi_probability_tsv, 
    generate_af_evolution_tsv,
)


def validate_file(filepath, file_type=None):
    """
    Validation function for argparse file arguments.
    """
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"File '{filepath}' does not exist.")
    if not os.access(filepath, os.R_OK):
        raise argparse.ArgumentTypeError(f"File '{filepath}' is not readable.")
    if file_type == "vcf" and not filepath.endswith((".vcf", ".vcf.gz")):
        raise argparse.ArgumentTypeError(f"File '{filepath}' is not a VCF file.")
    if file_type == "bed" and not filepath.endswith(".bed"):
        raise argparse.ArgumentTypeError(f"File '{filepath}' is not a BED file.")

    return filepath


def parse_args():
    """
    Parse command line arguments for the MSI analysis script.
    """
    parser = argparse.ArgumentParser(
        description="Quantify MSI Analysis: Compare MS Bed Regions with Varlociraptor VCF"
    )
    parser.add_argument(
        "--bed",
        type=lambda x: validate_file(x, "bed"),
        required=True,
        help="Path to the MSI bed file containing regions of interest.",
    )
    parser.add_argument(
        "--vcf",
        type=lambda x: validate_file(x, "vcf"),
        required=True,
        help="Path to the Varlociraptor VCF file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file path for the analysis results.",
    )
    parser.add_argument(
        "--debug-log",
        type=str,
        default="msi_analysis_debug.log",
        help="Output debug log file path (default: msi_analysis_debug.log)",
    )
    parser.add_argument(
        "--html-report",
        type=str,
        default="msi_report.html",
        help="Output HTML report file path (default: msi_report.html)",
    )
    parser.add_argument(
        "--quantification-output",
        type=str,
        default=None,
        help="Output file path for detailed quantification JSON data (optional)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for parallel processing",
    )
    parser.add_argument(
        "--unstable-threshold",
        type=float,
        default=0.5,
        help="P(≥1 MSI) threshold for calling region unstable (default: 0.5)"
    )
    parser.add_argument(
        "--msi-high-threshold", 
        type=float,
        default=3.5,
        help="MSI score threshold for MSI-High classification (default: 3.5 - MSIsensor standard)"
    )
    parser.add_argument(
        "--msi-probability-tsv",
        type=str,
        default=None,
        help="Output TSV file for MSI probability distribution (k, MSI score, probability)"
    )
    parser.add_argument(
        "--af-evolution-tsv", 
        type=str,
        default=None,
        help="Output TSV file for AF evolution summary across thresholds"
    )

    return parser.parse_args()


def main():
    """
    Main function to execute the MSI analysis.
    """

    args = parse_args()

    initialize_debug_log(args.debug_log, args)

    print("Loading VCF variants...")
    variants, filter_stats = load_vcf_variants(args.vcf)

    print("Loading BED regions...")
    regions, total_ms_regions = load_bed_regions(args.bed)

    print(f"Total MS regions from BED: {total_ms_regions:,}")


    print(f"Running intersection analysis with {args.threads} threads...")
    results, total_regions_loaded, unprocessed_count, merged_count = (
        bin_based_threaded_intersection(variants, regions, args.threads)
    )

    print("Preparing variants for DP analysis...")
    dp_result = prepare_variants_for_dp(results, args.vcf)

    
    if "error" not in dp_result and dp_result.get("ready_for_dp"):        
        print("Running AF evolution analysis...")
    
        af_evolution_results = run_af_evolution_analysis(
            dp_result,
            total_ms_regions,
            unstable_threshold=args.unstable_threshold,
            msi_high_threshold=args.msi_high_threshold,
        )

    else:
        print("Cannot run regional DP - data preparation failed")

    #TODO: Implement debug AF filtering discrepancy
    # debug_af_filtering_discrepancy(dp_result, total_ms_regions)

    output_data = {
        "analysis_info": {
            "timestamp": datetime.now().isoformat(),
            "vcf_file": args.vcf,
            "bed_file": args.bed,
            "threads_used": args.threads,
        },
        # "region_results": results,
        "region_results": dp_result,
        "af_evolution_results": af_evolution_results,
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=2, default=str)

    print(f"MSI analysis complete! Results saved to {args.output}")

    #TODO:ACTIVATE DEBUG LOGGING
    # write_complete_debug_log(
    #     args.debug_log,
    #     variants,
    #     regions,
    #     results,
    #     total_regions_loaded,
    #     args,
    #     filter_stats,
    #     unprocessed_count,
    #     merged_count,
    #     {"regional_analysis": regional_results,
    #     "af_evolution_results": af_evolution_results,
    #     },
    # )

    msi_data = create_msi_quantification(results)

    if args.quantification_output:
        with open(args.quantification_output, "w") as f:
            json.dump(msi_data, f, indent=2)
    print(f"Quantification data saved to: {args.quantification_output}")

    generate_msi_html_report(af_evolution_results, msi_data, args.html_report, args.msi_high_threshold)

    if args.msi_probability_tsv:
        generate_msi_probability_tsv(af_evolution_results, args.msi_probability_tsv)

    if args.af_evolution_tsv:
        generate_af_evolution_tsv(af_evolution_results, args.af_evolution_tsv)

    print(f"Debug log: {args.debug_log}")
    print(f"HTML report: {args.html_report}")


if __name__ == "__main__":
    main()
