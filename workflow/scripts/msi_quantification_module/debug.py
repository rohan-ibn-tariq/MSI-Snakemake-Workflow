"""
MSI Quantification Module - Debug Functions

Debug logging, data validation, and detailed analysis reporting.
Provides comprehensive debugging capabilities for MSI analysis workflow.
"""

from msi_quantification_module.core import create_msi_quantification
from datetime import datetime

def initialize_debug_log(debug_file_path, args):
    """Initialize debug log with header information"""
    with open(debug_file_path, "w") as debug_file:
        debug_file.write("=" * 80 + "\n")
        debug_file.write("MSI MICROSATELLITE INSTABILITY ANALYSIS - DEBUG LOG\n")
        debug_file.write("=" * 80 + "\n")
        debug_file.write(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        debug_file.write(f"VCF file: {args.vcf}\n")
        debug_file.write(f"BED file: {args.bed}\n")
        debug_file.write(f"Output JSON: {args.output}\n")
        debug_file.write(f"Debug log: {args.debug_log}\n")
        debug_file.write(f"HTML report: {args.html_report}\n")
        debug_file.write(f"Threads: {args.threads}\n")
        debug_file.write("=" * 80 + "\n\n")


def debug_data_loading_phase(debug_file, vcf_stats, total_regions_loaded):
    """Debug section for data loading phase"""
    debug_file.write("=" * 80 + "\n")
    debug_file.write("1. DATA LOADING DEBUG\n")
    debug_file.write("=" * 80 + "\n\n")

    debug_file.write("VCF PROCESSING FLOW:\n")
    debug_file.write(
        f"  Total VCF records in file: {vcf_stats['filter_stats']['total_records']:,}\n"
    )
    debug_file.write(
        f"  Records passed SNV filter: {vcf_stats['filter_stats']['records_passed']:,}\n"
    )
    debug_file.write(
        f"  Final indel variants processed: {vcf_stats['filter_stats']['variants_processed']:,}\n"
    )
    debug_file.write(
        f"  Difference: {vcf_stats['filter_stats']['records_passed'] - vcf_stats['filter_stats']['variants_processed']} (multi-allelic SNVs + problematic variants)\n\n"
    )

    debug_file.write("VCF VARIANT LOADING:\n")
    debug_file.write(f"  Total variant entries: {vcf_stats['total_count']:,}\n")
    debug_file.write(f"  Unique variant IDs: {vcf_stats['unique_count']:,}\n")
    debug_file.write(f"  Duplicates: {vcf_stats['duplicates']}\n\n")

    debug_file.write("EXCLUDED VARIANTS:\n")
    excluded_reasons = ["all_SNVs", "symbolic", "breakend", "spanning_deletion", "SNV"]
    for reason in excluded_reasons:
        if (
            reason in vcf_stats["filter_stats"]
            and vcf_stats["filter_stats"][reason] > 0
        ):
            debug_file.write(f"  {reason}: {vcf_stats['filter_stats'][reason]:,}\n")
    debug_file.write("\n")

    debug_file.write("BED REGION LOADING:\n")
    debug_file.write(f"  Total regions loaded: {total_regions_loaded:,}\n\n")


def debug_intersection_phase(debug_file, intersection_stats):
    """Debug section for intersection phase"""
    debug_file.write("=" * 80 + "\n")
    debug_file.write("2. INTERSECTION PROCESSING DEBUG\n")
    debug_file.write("=" * 80 + "\n\n")

    debug_file.write("BIN-BASED INTERSECTION SETUP:\n")
    debug_file.write(
        f"  Multi-threaded intersection with {intersection_stats['threads']} threads\n"
    )
    debug_file.write(f"  Bin tasks created: {intersection_stats['tasks']}\n")
    debug_file.write(
        f"  Total regions loaded: {intersection_stats['total_regions']:,}\n\n"
    )

    debug_file.write("INTERSECTION RESULTS:\n")
    debug_file.write(
        f"  Regions processed: {intersection_stats['regions_processed']:,}\n"
    )
    debug_file.write(
        f"  Regions with variants: {intersection_stats['regions_with_variants']:,}\n\n"
    )


def debug_numerical_reconciliation(debug_file, reconciliation_data):
    """Debug section for numerical reconciliation"""
    debug_file.write("=" * 80 + "\n")
    debug_file.write("3. NUMERICAL RECONCILIATION - DEBUG LOG\n")
    debug_file.write("=" * 80 + "\n\n")

    debug_file.write("ANALYSIS COUNT VERIFICATION:\n")
    debug_file.write(
        f"  Regions with variants: {reconciliation_data['regions_with_variants']:,}\n"
    )
    debug_file.write(
        f"  Total variant analyses: {reconciliation_data['total_analyses']:,}\n"
    )
    debug_file.write(f"  Difference: +{reconciliation_data['difference']}\n\n")

    debug_file.write("MULTI-VARIANT REGIONS EXPLANATION:\n")
    debug_file.write(
        f"  Regions with >1 variant: {reconciliation_data['multi_region_count']}\n"
    )
    debug_file.write("  BREAKDOWN:\n")
    for region in reconciliation_data["multi_regions"]:
        debug_file.write(
            f"    {region['region_id']}: {region['variant_count']} variants (+{region['extra_analyses']} extra)\n"
        )

    debug_file.write(
        f"\n  VERIFICATION: {reconciliation_data['regions_with_variants']} base + {reconciliation_data['total_extra']} extra = {reconciliation_data['total_analyses']} total\n\n"
    )


def debug_msi_analysis_results(debug_file, msi_analysis_data):
    """Debug section for MSI analysis results"""
    debug_file.write("=" * 80 + "\n")
    debug_file.write("4. MSI ANALYSIS RESULTS\n")
    debug_file.write("=" * 80 + "\n\n")

    debug_file.write("MSI REGION CLASSIFICATION:\n")
    debug_file.write(
        f"  Total microsatellite regions: {msi_analysis_data['total_regions']:,}\n"
    )
    debug_file.write(
        f"  MSI Unstable regions: {msi_analysis_data['msi_regions']:,} ({msi_analysis_data['msi_rate']:.2f}%)\n"
    )
    debug_file.write(
        f"  Uncertain regions (N/A indels): {msi_analysis_data['uncertain_regions']:,}\n"
    )
    debug_file.write(
        f"  MSI Stable regions: {msi_analysis_data['stable_regions']:,}\n\n"
    )

    debug_file.write("PERFECT MSI VARIANT BREAKDOWN:\n")
    debug_file.write(
        f"  Total perfect MSI variants: {msi_analysis_data['perfect_variants']:,}\n"
    )
    debug_file.write(
        f"  Injected MSI variants (MSI_*): {msi_analysis_data['injected_count']:,}\n"
    )
    debug_file.write(
        f"  Mason-created MSI variants: {msi_analysis_data['mason_count']:,}\n"
    )
    debug_file.write(f"  N/A variants: {msi_analysis_data['na_variants']:,}\n\n")

    debug_file.write("INDEL TYPE ANALYSIS:\n")
    debug_file.write(f"  Perfect MSI insertions: {msi_analysis_data['insertions']:,}\n")
    debug_file.write(f"  Perfect MSI deletions: {msi_analysis_data['deletions']:,}\n")
    debug_file.write(
        f"  Insertion/Deletion ratio: {msi_analysis_data['indel_ratio']:.2f}\n\n"
    )

    debug_file.write("REGION vs VARIANT RECONCILIATION:\n")
    debug_file.write(
        f"  Perfect MSI variants: {msi_analysis_data['perfect_variants']:,}\n"
    )
    debug_file.write(f"  MSI unstable regions: {msi_analysis_data['msi_regions']:,}\n")
    debug_file.write(f"  Difference: {msi_analysis_data['total_extra_perfect']}\n")

    if msi_analysis_data["regions_with_multiple_perfect"]:
        debug_file.write(
            f"  Regions with multiple perfect variants: {len(msi_analysis_data['regions_with_multiple_perfect'])}\n"
        )
        for region in msi_analysis_data["regions_with_multiple_perfect"]:
            debug_file.write(
                f"    {region['region_id']}: {region['perfect_count']} perfect variants\n"
            )
            for var_id in region["variant_ids"]:
                debug_file.write(f"      - {var_id}\n")
    debug_file.write("\n")


def debug_quantification_summary(debug_file, quantification_data):
    """Debug section for quantification summary"""
    debug_file.write("=" * 80 + "\n")
    debug_file.write("5. MSI QUANTIFICATION SUMMARY\n")
    debug_file.write("=" * 80 + "\n\n")

    debug_file.write("ALLELE FREQUENCY ANALYSIS (Perfect MSI only):\n")
    debug_file.write(
        f"  Variants with AF data: {quantification_data['af_with_data']:,}\n"
    )
    debug_file.write(f"  Variants with zero AF: {quantification_data['af_zero']:,}\n")
    debug_file.write(f"  Variants with N/A AF: {quantification_data['af_na']:,}\n")
    debug_file.write(
        f"  Mean AF (non-zero, non-N/A): {quantification_data['af_mean']:.3f}\n\n"
    )

    debug_file.write("MOTIF TYPE MSI RATES:\n")
    for motif_type, data in sorted(quantification_data["motif_breakdown"].items()):
        msi_rate = (data["msi"] / data["total"] * 100) if data["total"] > 0 else 0
        debug_file.write(
            f"  {motif_type}: {data['msi']:,}/{data['total']:,} ({msi_rate:.1f}% MSI rate)\n"
        )

    debug_file.write("\nDETAILED MOTIF BREAKDOWN:\n")
    for motif_type, data in sorted(quantification_data["motif_breakdown"].items()):
        debug_file.write(f"  {motif_type} ({data['total']:,} total):\n")
        debug_file.write(f"    MSI Unstable: {data['msi']:,}\n")
        debug_file.write(
            f"    Uncertain (N/A indels): {data['non_msi_with_indels']:,}\n"
        )
        debug_file.write(f"    MSI Stable: {data['stable']:,}\n")
    debug_file.write("\n")


def print_msi_summary(msi_data):
    print("\n" + "=" * 60)
    print("MSI QUANTIFICATION SUMMARY")
    print("=" * 60)

    scope = msi_data["data_scope"]
    print(f"\nDATA SCOPE:")
    print(f"  Total microsatellite regions: {scope['total_ms_regions']:,}")
    print(
        f"  MSI regions (perfect indels): {scope['msi_regions_with_perfect_indels']:,}"
    )
    print(f"  Non-MSI regions with indels: {scope['non_msi_regions_with_indels']:,}")
    print(f"  Stable regions: {scope['stable_regions']:,}")
    print(f"  Total variants found: {scope['total_variants_found']:,}")
    print(f"  Perfect MSI variants: {scope['perfect_msi_variants_analyzed']:,}")
    print(f"  N/A variants: {scope['na_variants_excluded']:,}")

    overview = msi_data["overview"]
    print(f"\nMSI ANALYSIS:")
    print(f"  MSI regions: {overview['msi_regions']:,} ({overview['msi_rate']:.2f}%)")
    print(f"  Non-MSI with indels: {overview['non_msi_regions_with_indels']:,}")
    print(f"  Stable regions: {overview['stable_regions']:,}")
    print(f"  Total: {overview['total_regions']:,}")

    indels = msi_data["indels"]
    print(f"\nPERFECT MSI INDELS:")
    print(f"  Insertions: {indels['insertion_count']:,}")
    print(f"  Deletions: {indels['deletion_count']:,}")
    print(
        f"  Total perfect MSI indels: {indels['insertion_count'] + indels['deletion_count']:,}"
    )
    print(
        f"  Insertion/Deletion ratio: {indels['insertion_count']/indels['deletion_count']:.2f}"
        if indels["deletion_count"] > 0
        else "  No deletions found"
    )

    af = msi_data["allele_frequencies"]
    print(f"\nALLELE FREQUENCY (Perfect MSI only):")
    print(f"  Variants with AF data: {af['count']:,}")
    print(f"  Variants with zero AF: {af['count_with_zero_af']:,}")
    print(f"  Variants with N/A AF: {af['count_with_na_af']:,}")
    print(f"  Mean AF (non-zero, non-N/A): {af['mean']:.3f}")

    print(f"\nMOTIF TYPE BREAKDOWN:")
    motif_data = msi_data["motif_breakdown"]
    for motif_type, data in sorted(motif_data.items()):
        msi_rate = (data["msi"] / data["total"]) * 100 if data["total"] > 0 else 0
        print(f"  {motif_type}: {data['total']:,} total")
        print(f"    MSI: {data['msi']:,} ({msi_rate:.1f}%)")
        print(f"    Non-MSI with indels: {data['non_msi_with_indels']:,}")
        print(f"    Stable: {data['stable']:,}")

    print(f"\nNote: {scope['analysis_note']}")
    print("=" * 60)


def debug_dp_analysis_results(debug_file, dp_data):
    """Debug section for DP-based MSI analysis results."""
    
    debug_file.write("\n" + "="*80 + "\n")
    debug_file.write("6. DP-BASED MSI ANALYSIS RESULTS\n")
    debug_file.write("="*80 + "\n")
    
    regional = dp_data.get("regional_analysis", {})
    if regional:
        debug_file.write("\nREGIONAL MSI ANALYSIS (DP Method):\n")
        debug_file.write(f"  MSI Score: {regional.get('msi_score', 'N/A')}% ± {regional.get('msi_uncertainty', 'N/A')}%\n")
        debug_file.write(f"  MSI Status: {regional.get('msi_status', 'N/A')}\n")
        debug_file.write(f"  Unstable Regions: {regional.get('unstable_regions', 0):,}\n")
        debug_file.write(f"  Total MS Regions: {regional.get('total_regions', 0):,}\n")
        debug_file.write(f"  Regions with Variants: {regional.get('regions_with_variants', 0):,}\n")
        debug_file.write(f"  Regions without Variants: {regional.get('regions_without_variants', 0):,}\n")
        debug_file.write(f"  Expected Unstable Regions: {regional.get('expected_unstable_regions', 'N/A')} ± {regional.get('expected_unstable_uncertainty', 'N/A')}\n")
        debug_file.write(f"  Expected MSI Variants: {regional.get('expected_msi_variants', 'N/A')} ± {regional.get('expected_variants_uncertainty', 'N/A')}\n")


        total_regions_value = regional.get('total_regions', 0)
        regions_with_variants = regional.get('regions_with_variants', 0)
        
        if total_regions_value > 0:
            coverage = (regions_with_variants / total_regions_value) * 100
            debug_file.write(f"  Analysis Coverage: {coverage:.1f}%\n")
        else:
            debug_file.write(f"  Analysis Coverage: ERROR - total_regions = {total_regions_value}\n")
        debug_file.write("\n")
    
    af_evolution = dp_data.get("af_evolution", {})
    if af_evolution:
        debug_file.write("\nAF EVOLUTION TIMELINE:\n")
        for af_key in sorted(af_evolution.keys()):
            af_data = af_evolution[af_key]
            debug_file.write(f"  {af_key}: {af_data.get('msi_score', 'N/A')}%")
            if 'msi_uncertainty' in af_data:
                debug_file.write(f" ± {af_data['msi_uncertainty']}%")
            debug_file.write("\n")


def write_complete_debug_log(
    debug_file_path,
    variants,
    regions,
    results,
    total_regions_loaded,
    args,
    filter_stats,
    unprocessed_count,
    merged_count,
    dp_results=None,
):
    """Master debug function - handles header and all debug sections"""

    initialize_debug_log(debug_file_path, args)

    ### Collect VCF Stats ###
    total_variants = sum(
        len(bin_variants)
        for chrom_bins in variants.values()
        for bin_variants in chrom_bins.values()
    )
    unique_variant_ids = set()
    for chrom_bins in variants.values():
        for bin_variants in chrom_bins.values():
            for variant_id, variant_data in bin_variants:
                unique_variant_ids.add(variant_id)

    vcf_stats = {
        "total_count": total_variants,
        "unique_count": len(unique_variant_ids),
        "duplicates": total_variants - len(unique_variant_ids),
        "filter_stats": filter_stats,
    }
    #### END VCF STATS #####

    ### Collect Intersection Stats ###
    tasks_count = 0
    for chrom in regions:
        for bin_id in regions[chrom]:
            if chrom in variants and bin_id in variants[chrom]:
                tasks_count += 1

    regions_with_variants = sum(1 for r in results if r.get("has_variants", False))
    total_analyses = sum(len(r.get("variants", [])) for r in results)

    variant_analysis_count = {}
    for r in results:
        for v in r.get("variants", []):
            variant_id = v["variant_id"]
            if variant_id not in variant_analysis_count:
                variant_analysis_count[variant_id] = 0
            variant_analysis_count[variant_id] += 1

    unique_variants_in_results = set(variant_analysis_count.keys())

    multi_variant_regions = []
    for r in results:
        variant_count = len(r.get("variants", []))
        if variant_count > 1:
            multi_variant_regions.append(
                {
                    "region_id": r.get("region_id"),
                    "variant_count": variant_count,
                    "extra_analyses": variant_count - 1,
                }
            )

    total_extra_analyses = sum(
        region["extra_analyses"] for region in multi_variant_regions
    )

    intersection_stats = {
        "threads": args.threads,
        "tasks": tasks_count,
        "total_regions": total_regions_loaded,
        "unprocessed_regions": unprocessed_count,
        "merged_results": merged_count,
        "regions_processed": len(results),
        "regions_with_variants": regions_with_variants,
        "total_analyses": total_analyses,
        "unique_variants": len(unique_variants_in_results),
        "avg_analyses_per_variant": (
            total_analyses / len(unique_variants_in_results)
            if unique_variants_in_results
            else 0
        ),
        "multi_variant_regions": multi_variant_regions,
        "total_extra_analyses": total_extra_analyses,
    }
    #### END Intersection STATS ######

    #### NUMERICAL RECONCILIATION ####
    reconciliation_data = {
        "regions_with_variants": regions_with_variants,
        "total_analyses": total_analyses,
        "difference": total_analyses - regions_with_variants,
        "multi_region_count": len(multi_variant_regions),
        "multi_regions": multi_variant_regions,
        "total_extra": total_extra_analyses,
    }
    #### END NUMERICAL RECONCILIATION ####

    #### MSI ANALYSIS DATA ####
    msi_regions = sum(1 for r in results if r.get("num_perfect_repeats", 0) > 0)
    uncertain_regions = sum(1 for r in results if r.get("num_na_repeats", 0) > 0)
    stable_regions = len(results) - msi_regions - uncertain_regions
    perfect_variants = sum(r.get("num_perfect_repeats", 0) for r in results)
    na_variants = sum(r.get("num_na_repeats", 0) for r in results)

    injected_count = 0
    mason_count = 0
    insertions = 0
    deletions = 0

    for r in results:
        for v in r.get("variants", []):
            if v.get("repeat_status") == "perfect":
                vcf_id = v.get("vcf_id", "")
                if vcf_id and vcf_id.startswith("MSI_"):
                    injected_count += 1
                else:
                    mason_count += 1

                if v["variant_type"] == "insertion":
                    insertions += 1
                else:
                    deletions += 1

    regions_with_multiple_perfect = []
    for r in results:
        perfect_count = r.get("num_perfect_repeats", 0)
        if perfect_count > 1:
            perfect_variant_ids = []
            for v in r.get("variants", []):
                if v.get("repeat_status") == "perfect":
                    perfect_variant_ids.append(v["variant_id"])

            regions_with_multiple_perfect.append(
                {
                    "region_id": r.get("region_id"),
                    "perfect_count": perfect_count,
                    "variant_ids": perfect_variant_ids,
                }
            )

    total_extra_perfect = sum(
        region["perfect_count"] - 1 for region in regions_with_multiple_perfect
    )

    msi_analysis_data = {
        "total_regions": len(results),
        "msi_regions": msi_regions,
        "msi_rate": (msi_regions / len(results)) * 100,
        "uncertain_regions": uncertain_regions,
        "stable_regions": stable_regions,
        "perfect_variants": perfect_variants,
        "injected_count": injected_count,
        "mason_count": mason_count,
        "na_variants": na_variants,
        "insertions": insertions,
        "deletions": deletions,
        "indel_ratio": insertions / deletions if deletions > 0 else float("inf"),
        "regions_with_multiple_perfect": regions_with_multiple_perfect,
        "total_extra_perfect": total_extra_perfect,
    }
    #### END MSI ANALYSIS DATA ####

    #### QUANTIFICATION DATA ######
    temp_msi_data = create_msi_quantification(results)

    quantification_data = {
        "af_with_data": temp_msi_data["allele_frequencies"]["count"],
        "af_zero": temp_msi_data["allele_frequencies"]["count_with_zero_af"],
        "af_na": temp_msi_data["allele_frequencies"]["count_with_na_af"],
        "af_mean": temp_msi_data["allele_frequencies"]["mean"],
        "motif_breakdown": temp_msi_data["motif_breakdown"],
    }
    #### END QUANTIFICATION DATA ####

    with open(debug_file_path, "a") as debug_file:
        debug_data_loading_phase(debug_file, vcf_stats, total_regions_loaded)
        debug_intersection_phase(debug_file, intersection_stats)
        debug_numerical_reconciliation(debug_file, reconciliation_data)
        debug_msi_analysis_results(debug_file, msi_analysis_data)
        debug_quantification_summary(debug_file, quantification_data)
        if dp_results:
            debug_dp_analysis_results(debug_file, dp_results)

