"""
Quantify MSI Analysis: 
    Compare MS Bed Regions with Varlociraptor VCF
"""

#!/usr/bin/env python3

import argparse
import json
import multiprocessing
import os
from collections import defaultdict
from datetime import datetime

import altair as alt
import polars as pl

import pysam
from binning import assign_bin

alt.data_transformers.disable_max_rows()


#################################################DEBUG  FUNCTIONS###############################
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


#######################################################################################################

###MAIN END
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
###MAIN END

####CORE START
def normalize_chromosome(chrom):
    """
    Normalize chromosome names to consistent format
    """
    if chrom.startswith("chr"):
        return chrom
    return "chr" + chrom


def calculate_dynamic_svlen_with_anchor(ref, alt):
    """Calculate SVLEN using dynamic anchor detection (matches analysis logic)"""

    anchor_len = 0
    min_len = min(len(ref), len(alt))

    for i in range(min_len):
        if ref[i].upper() == alt[i].upper():
            anchor_len += 1
        else:
            break

    if len(alt) > len(ref):
        changed_seq = alt[anchor_len:]
        svlen = len(changed_seq)
    else:
        changed_seq = ref[anchor_len:]
        svlen = -len(changed_seq)

    return svlen


def load_bed_regions(bed_file):
    """
    Load all microsatellite regions from UCSC BED file.
    Parse motif information and calculate buffer sizes.
    """
    regions = defaultdict(lambda: defaultdict(list))

    print(f"[MSI-ANALYSIS INFO] Loading BED regions from: {bed_file}")

    with open(bed_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                print(
                    f"[MSI-ANALYSIS WARNING] Skipping line {line_num}: insufficient fields"
                )
                continue

            try:
                bin_id = int(fields[0])
                chrom = fields[1]
                start = int(fields[2])
                end = int(fields[3])
                motif_string = fields[4]

                normalized_chrom = normalize_chromosome(chrom)

                if "x" not in motif_string:
                    print(
                        f"[MSI-ANALYSIS WARNING] Invalid motif format at "
                        f"line {line_num}: {motif_string}"
                    )
                    continue

                count_str, motif = motif_string.split("x", 1)
                repeat_count = int(count_str)
                motif_length = len(motif)
                buffer_size = motif_length

                motif_types = {
                    1: "mononucleotide",
                    2: "dinucleotide",
                    3: "trinucleotide",
                    4: "tetranucleotide",
                    5: "pentanucleotide",
                    6: "hexanucleotide",
                }
                motif_type = motif_types.get(motif_length, f"{motif_length}mer")

                region_id = f"{normalized_chrom}:{start}-{end}"

                region_data = {
                    "chrom": normalized_chrom,
                    "start": start,
                    "end": end,
                    "motif_string": motif_string,
                    "repeat_count": repeat_count,
                    "motif": motif,
                    "motif_length": motif_length,
                    "motif_type": motif_type,
                    "buffer_size": buffer_size,
                    "region_id": region_id,
                }

                regions[normalized_chrom][bin_id].append((region_id, region_data))

            except (ValueError, IndexError) as e:
                print(f"[MSI-ANALYSIS WARNING] Error parsing line {line_num}: {e}")
                continue

    if not regions:
        raise ValueError(f"[MSI-ANALYSIS ERROR] No valid regions found in {bed_file}")

    total_regions = sum(
        len(bin_regions)
        for chrom_bins in regions.values()
        for bin_regions in chrom_bins.values()
    )
    print(f"[MSI-ANALYSIS INFO] Loaded {total_regions} microsatellite regions")
    return regions


def is_perfect_repeat(sequence, motif):
    """Check if sequence is a perfect repeat of the motif (complete motifs only)"""
    if not sequence or not motif:
        return False

    remainder = len(sequence) % len(motif)
    if remainder > 0:
        return False

    repeat_count = len(sequence) // len(motif)
    expected = motif * repeat_count

    return sequence.upper() == expected.upper()


def check_variant_overlaps_region(variant_data, region_data):
    """Check if variant overlaps with MS region (inclusive boundaries)"""
    variant_start = variant_data["pos"] - 1

    if len(variant_data["alt"]) > len(variant_data["ref"]):
        variant_end = variant_start + 1
    else:
        variant_end = variant_start + 1

    region_start = region_data["start"]
    region_end = region_data["end"]

    overlaps = not (variant_end <= region_start or variant_start >= region_end)

    return overlaps


def analyze_variant_in_region(variant_data, region_data):
    """Analyze a variant within a microsatellite region"""

    if not check_variant_overlaps_region(variant_data, region_data):
        return None

    ref_seq = variant_data["ref"]
    alt_seq = variant_data["alt"]

    anchor_len = 0
    min_len = min(len(ref_seq), len(alt_seq))

    for i in range(min_len):
        if ref_seq[i].upper() == alt_seq[i].upper():
            anchor_len += 1
        else:
            break

    if anchor_len == 0:
        print(
            f"[WARNING] No anchor base found for variant {variant_data['variant_id']}"
        )
        return None

    if variant_data["svlen"] > 0:
        changed_seq = alt_seq[anchor_len:]
    else:
        changed_seq = ref_seq[anchor_len:]

    if not changed_seq:
        print(f"[WARNING] No changed sequence for variant {variant_data['variant_id']}")
        return None

    motif = region_data["motif"]
    is_perfect = is_perfect_repeat(changed_seq, motif)
    repeat_status = "perfect" if is_perfect else "N/A"


    return {
        "variant_id": variant_data["variant_id"],
        "vcf_id": variant_data.get("vcf_id"),
        "region_id": region_data["region_id"],
        "motif": motif,
        "motif_type": region_data["motif_type"],
        "variant_type": variant_data["variant_type"],
        "svlen": variant_data["svlen"],
        "changed_sequence": changed_seq,
        "repeat_status": repeat_status,
        "af_mean": variant_data["af_mean"],
        "af_max": variant_data["af_max"],
    }


def load_vcf_variants(vcf_file):
    """
    Load all indel variants from VCF file using pysam.
    Handle multi-allelic sites - creates separate variant records for each alt.
    """
    variants = defaultdict(lambda: defaultdict(list))
    filter_stats = defaultdict(int)

    print(f"[MSI-ANALYSIS INFO] Loading VCF variants from: {vcf_file}")

    try:
        vcf = pysam.VariantFile(vcf_file)

        for record in vcf:
            filter_stats["total_records"] += 1

            if all(len(alt) == len(record.ref) for alt in record.alts):
                filter_stats["all_SNVs"] += 1
                continue

            filter_stats["records_passed"] += 1

            prob_present_info = record.info.get("PROB_PRESENT", None)
            prob_absent_info = record.info.get("PROB_ABSENT", None)
            prob_artifact_info = record.info.get("PROB_ARTIFACT", None)

            if prob_present_info is None:
                prob_present_list = [None]
            elif isinstance(prob_present_info, (list, tuple)):
                prob_present_list = list(prob_present_info)
            else:
                prob_present_list = [prob_present_info]

            if prob_absent_info is None:
                prob_absent_list = [None]
            elif isinstance(prob_absent_info, (list, tuple)):
                prob_absent_list = list(prob_absent_info)
            else:
                prob_absent_list = [prob_absent_info]

            if prob_artifact_info is None:
                prob_artifact_list = [None]
            elif isinstance(prob_artifact_info, (list, tuple)):
                prob_artifact_list = list(prob_artifact_info)
            else:
                prob_artifact_list = [prob_artifact_info]

            svlen_info = record.info.get("SVLEN")
            svlen_list = (
                list(svlen_info)
                if isinstance(svlen_info, (list, tuple))
                else [svlen_info] if svlen_info is not None else None
            )

            all_sample_afs = {}
            for sample_name, sample_data in record.samples.items():
                if "AF" in sample_data:
                    all_sample_afs[sample_name] = sample_data["AF"]

            for alt_idx, alt in enumerate(record.alts):
                if len(alt) == len(record.ref):
                    filter_stats["SNV"] += 1
                    continue

                if "<" in alt or ">" in alt or "<" in record.ref or ">" in record.ref:
                    filter_stats["symbolic"] += 1
                    continue

                if "[" in alt or "]" in alt:
                    filter_stats["breakend"] += 1
                    continue

                if alt == "*":
                    filter_stats["spanning_deletion"] += 1
                    continue

                filter_stats["variants_processed"] += 1

                prob_present = (
                    prob_present_list[alt_idx]
                    if alt_idx < len(prob_present_list)
                    else None
                )
                prob_absent = (
                    prob_absent_list[alt_idx]
                    if alt_idx < len(prob_absent_list)
                    else None
                )
                prob_artifact = (
                    prob_artifact_list[alt_idx]
                    if alt_idx < len(prob_artifact_list)
                    else None
                )

                dynamic_svlen = calculate_dynamic_svlen_with_anchor(record.ref, alt)

                if (
                    svlen_list is not None
                    and alt_idx < len(svlen_list)
                    and svlen_list[alt_idx] is not None
                ):
                    svlen = svlen_list[alt_idx]
                else:
                    svlen = dynamic_svlen

                variant_id = f"{record.chrom}:{record.pos}:{alt_idx}"

                all_sample_afs_for_alt = {}
                for sample_name, af_data in all_sample_afs.items():
                    af_value = None
                    if isinstance(af_data, (list, tuple)):
                        if alt_idx < len(af_data) and af_data[alt_idx] is not None:
                            af_value = float(af_data[alt_idx])
                    else:
                        if alt_idx == 0 and af_data is not None:
                            af_value = float(af_data)

                    if af_value is not None:
                        all_sample_afs_for_alt[sample_name] = af_value

                non_zero_afs = {
                    sample: af
                    for sample, af in all_sample_afs_for_alt.items()
                    if af > 0.0
                }

                if non_zero_afs:
                    af_mean = sum(non_zero_afs.values()) / len(non_zero_afs)
                    af_max = max(non_zero_afs.values())
                    af_min = min(non_zero_afs.values())
                elif all_sample_afs_for_alt:
                    af_mean = 0.0
                    af_max = 0.0
                    af_min = 0.0
                else:
                    af_mean = "N/A"
                    af_max = "N/A"
                    af_min = "N/A"

                num_samples_with_variant = len(non_zero_afs)
                num_samples_total = len(all_sample_afs_for_alt)
                variant_penetrance = (
                    num_samples_with_variant / num_samples_total
                    if num_samples_total > 0
                    else 0.0
                )

                normalized_chrom = normalize_chromosome(record.chrom)
                variant_type = "insertion" if svlen > 0 else "deletion"

                variant_data = {
                    "chrom": normalized_chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": alt,
                    "vcf_id": record.id if record.id and record.id != "." else None,
                    "svlen": svlen,
                    "variant_type": variant_type,
                    "prob_present": (
                        10 ** (-float(prob_present) / 10)
                        if prob_present is not None and prob_present != float("inf")
                        else 0.0
                    ),
                    "prob_absent": (
                        10 ** (-float(prob_absent) / 10)
                        if prob_absent is not None and prob_absent != float("inf")
                        else 0.0
                    ),
                    "prob_artifact": (
                        10 ** (-float(prob_artifact) / 10)
                        if prob_artifact is not None and prob_artifact != float("inf")
                        else 0.0
                    ),
                    "sample_afs": all_sample_afs_for_alt,
                    "af_mean": af_mean,
                    "af_max": af_max,
                    "af_min": af_min,
                    "num_samples_with_variant": num_samples_with_variant,
                    "num_samples_total": num_samples_total,
                    "variant_penetrance": variant_penetrance,
                    "variant_id": variant_id,
                    "alt_index": alt_idx,
                }

                variant_start = record.pos - 1
                variant_end = variant_start + max(len(record.ref), len(alt))

                bin_id = assign_bin(variant_start, variant_end)
                variants[normalized_chrom][bin_id].append((variant_id, variant_data))

        vcf.close()

    except Exception as e:
        raise ValueError(f"[MSI-ANALYSIS ERROR] Error loading VCF: {e}")

    if not variants:
        print("[MSI-ANALYSIS WARNING] No indel variants found in VCF")
    else:
        total_count = sum(
            len(bin_variants)
            for chrom_bins in variants.values()
            for bin_variants in chrom_bins.values()
        )
        print(
            f"[MSI-ANALYSIS INFO] Loaded {total_count} indel variants with bins assigned"
        )
        print("[MSI-ANALYSIS INFO] Variant filtering results:")
        for reason, count in filter_stats.items():
            print(f"  {reason}: {count}")

    return variants, filter_stats


def create_bin_tasks(regions, variants):
    """Create list of tasks for parallel processing"""
    tasks = []
    total_regions_loaded = 0
    unprocessed_regions = []

    for chrom in regions:
        for bin_id in regions[chrom]:
            regions_in_bin = regions[chrom][bin_id]
            total_regions_loaded += len(regions_in_bin)

            if chrom in variants and bin_id in variants[chrom]:
                variants_in_bin = variants[chrom][bin_id]
                tasks.append((chrom, bin_id, regions_in_bin, variants_in_bin))
            else:
                unprocessed_regions.extend(regions_in_bin)

    print(f"[MSI-ANALYSIS INFO] Total regions loaded: {total_regions_loaded:,}")
    print(f"[MSI-ANALYSIS INFO] Created {len(tasks)} bin tasks for parallel processing")
    print(f"[MSI-ANALYSIS INFO] Unprocessed regions: {len(unprocessed_regions):,}")

    return tasks, total_regions_loaded, unprocessed_regions


def process_single_bin(chrom, bin_id, regions_in_bin, variants_in_bin):
    """Process one bin independently - core MSI intersection logic"""

    bin_results = []

    for region_id, region_data in regions_in_bin:
        variants_in_region = []

        for variant_id, variant_data in variants_in_bin:
            if "pos" not in variant_data:
                print("ERROR: variant_data missing 'pos' key!")
                print(f"ERROR: variant_data = {variant_data}")
                continue

            if variant_data["chrom"] != region_data["chrom"]:
                continue

            analysis = analyze_variant_in_region(variant_data, region_data)
            if analysis is not None:
                variants_in_region.append(analysis)

        region_summary = {
            "region_id": region_id,
            "chrom": region_data["chrom"],
            "start": region_data["start"],
            "end": region_data["end"],
            "motif": region_data["motif"],
            "motif_type": region_data["motif_type"],
            "has_variants": len(variants_in_region) > 0,
            "num_insertions": len(
                [v for v in variants_in_region if v["variant_type"] == "insertion"]
            ),
            "num_deletions": len(
                [v for v in variants_in_region if v["variant_type"] == "deletion"]
            ),
            "num_perfect_repeats": len(
                [v for v in variants_in_region if v["repeat_status"] == "perfect"]
            ),
            "num_na_repeats": len(
                [v for v in variants_in_region if v["repeat_status"] == "N/A"]
            ),
            "max_af": max([v["af_max"] for v in variants_in_region], default=0.0),
            "variants": variants_in_region,
        }

        bin_results.append(region_summary)

    return bin_results


def merge_bin_results(bin_results):
    """Merge results from all processed bins"""
    final_results = []

    for bin_result in bin_results:
        if bin_result:
            final_results.extend(bin_result)

    print(f"[MSI-ANALYSIS INFO] Merged {len(final_results)} region results")
    return final_results


def bin_based_threaded_intersection(variants, regions, n_threads):
    """
    Fast intersection using UCSC bins + threading
    """
    print(
        f"[MSI-ANALYSIS INFO] Starting multi-threaded intersection with {n_threads} threads"
    )

    tasks, total_regions_loaded, unprocessed_regions = create_bin_tasks(
        regions, variants
    )

    unprocessed_count = len(unprocessed_regions)

    with multiprocessing.Pool(n_threads) as pool:
        bin_results = pool.starmap(process_single_bin, tasks)

    final_results = merge_bin_results(bin_results)
    merged_count = len(final_results)

    for region_id, region_data in unprocessed_regions:
        empty_result = {
            "region_id": region_id,
            "chrom": region_data["chrom"],
            "start": region_data["start"],
            "end": region_data["end"],
            "motif": region_data["motif"],
            "motif_type": region_data["motif_type"],
            "has_variants": False,
            "num_insertions": 0,
            "num_deletions": 0,
            "num_perfect_repeats": 0,
            "num_na_repeats": 0,
            "max_af": 0.0,
            "variants": [],
        }
        final_results.append(empty_result)

    regions_unprocessed = total_regions_loaded - len(final_results)
    print("[MSI-ANALYSIS INFO] Completed intersection:")
    print(f"  Regions processed: {len(final_results):,}")
    print(f"  Regions unprocessed: {regions_unprocessed:,}")
    print(f"  Total regions: {total_regions_loaded:,}")
    return final_results, total_regions_loaded, unprocessed_count, merged_count


### CORE END

### MAIN START
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
        "--threads",
        type=int,
        default=4,
        help="Number of threads for parallel processing",
    )

    return parser.parse_args()
#MAIN END

#### ADDENDUM CORE ####
def create_msi_quantification(results):

    msi_regions = len([r for r in results if r.get("num_perfect_repeats", 0) > 0])
    non_msi_regions_with_indels = len(
        [
            r
            for r in results
            if r.get("has_variants") == True and r.get("num_perfect_repeats", 0) == 0
        ]
    )
    stable_regions = len([r for r in results if r.get("has_variants") == False])
    total_variants = sum(len(r.get("variants", [])) for r in results)
    perfect_variants = sum(r.get("num_perfect_repeats", 0) for r in results)
    na_variants = total_variants - perfect_variants

    insertion_sizes = []
    deletion_sizes = []
    insertion_units = []
    deletion_units = []
    allele_frequencies = []
    perfect_insertion_count = 0
    perfect_deletion_count = 0
    af_zero_count = 0
    af_na_count = 0

    for region in results:
        for variant in region.get("variants", []):
            if variant.get("repeat_status") == "perfect":
                if variant["variant_type"] == "insertion":
                    perfect_insertion_count += 1
                    insertion_sizes.append(variant["svlen"])
                    if len(variant["motif"]) > 0:
                        insertion_units.append(variant["svlen"] / len(variant["motif"]))
                else:
                    perfect_deletion_count += 1
                    deletion_sizes.append(abs(variant["svlen"]))
                    if len(variant["motif"]) > 0:
                        deletion_units.append(
                            abs(variant["svlen"]) / len(variant["motif"])
                        )

                if isinstance(variant["af_max"], (int, float)):
                    if variant["af_max"] > 0:
                        allele_frequencies.append(variant["af_max"])
                    else:
                        af_zero_count += 1
                else:
                    af_na_count += 1

    motif_summary = {}
    for region in results:
        motif_type = region.get("motif_type", "unknown")
        if motif_type not in motif_summary:
            motif_summary[motif_type] = {
                "total": 0,
                "msi": 0,
                "non_msi_with_indels": 0,
                "stable": 0,
            }

        motif_summary[motif_type]["total"] += 1
        if region.get("num_perfect_repeats", 0) > 0:
            motif_summary[motif_type]["msi"] += 1
        elif region.get("has_variants") == True:
            motif_summary[motif_type]["non_msi_with_indels"] += 1
        else:
            motif_summary[motif_type]["stable"] += 1

    return {
        "data_scope": {
            "total_ms_regions": len(results),
            "msi_regions_with_perfect_indels": msi_regions,
            "non_msi_regions_with_indels": non_msi_regions_with_indels,
            "stable_regions": stable_regions,
            "total_variants_found": total_variants,
            "perfect_msi_variants_analyzed": perfect_variants,
            "na_variants_excluded": na_variants,
            "analysis_note": "All histograms and AF analysis based on perfect MSI variants only",
        },
        "overview": {
            "total_regions": len(results),
            "msi_regions": msi_regions,
            "non_msi_regions_with_indels": non_msi_regions_with_indels,
            "stable_regions": stable_regions,
            "msi_rate": (msi_regions / len(results)) * 100,
        },
        "indels": {
            "insertion_count": perfect_insertion_count,
            "deletion_count": perfect_deletion_count,
            "insertion_sizes": insertion_sizes,
            "deletion_sizes": deletion_sizes,
            "insertion_units": insertion_units,
            "deletion_units": deletion_units,
        },
        "allele_frequencies": {
            "values": allele_frequencies,
            "count": len(allele_frequencies),
            "mean": (
                sum(allele_frequencies) / len(allele_frequencies)
                if allele_frequencies
                else 0
            ),
            "count_with_zero_af": af_zero_count,
            "count_with_na_af": af_na_count,
            "mean_note": "Mean calculated from non-zero, non-N/A values only",
            "note": "AF values from perfect MSI variants only",
        },
        "motif_breakdown": motif_summary,
    }
### END ADDENDUM CORE ###

###TO REMOVE:
def analyze_msi_results_post_quantification(results):
    print("[MSI-ANALYSIS] Creating quantification...")

    msi_data = create_msi_quantification(results)

    with open("msi_quantification.json", "w") as f:
        json.dump(msi_data, f, indent=2)

    print("[MSI-ANALYSIS] Quantification saved to: msi_quantification.json")

    return msi_data
#### TO REMOVE END ####

#########################################HTML REPORTING#########################################


def generate_msi_html_report(msi_data, output_path):
    """
    Generate MSI quantification HTML report.
    """

    overview_chart = create_overview_chart(msi_data)
    simple_indel_chart = create_simple_indel_chart(msi_data)
    indel_histogram = create_indel_histogram(msi_data)
    af_distribution = create_af_distribution_chart(msi_data)
    motif_breakdown_chart = create_motif_breakdown_chart(msi_data)

    scope = msi_data["data_scope"]
    overview = msi_data["overview"]
    indels = msi_data["indels"]
    af = msi_data["allele_frequencies"]

    timestamp = datetime.now().strftime("%B %d, %Y")

    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MSI Quantification Report</title>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <style>
        body {{
            font-family: "Times New Roman", Times, serif;
            background-color: #fefdfb;
            color: #1a1a1a;
            line-height: 1.6;
            max-width: 900px;
            margin: 0 auto;
            padding: 40px 20px;
        }}
        h1 {{
            text-align: center;
            font-size: 24px;
            font-weight: normal;
            border-bottom: 2px solid #1a1a1a;
            padding-bottom: 10px;
            margin-bottom: 30px;
        }}
        h2 {{
            font-size: 18px;
            font-weight: normal;
            margin-top: 40px;
            margin-bottom: 15px;
            border-bottom: 1px solid #666;
        }}
        .header-info {{
            text-align: center;
            font-style: italic;
            margin-bottom: 40px;
            color: #666;
        }}
        .metrics {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .metric {{
            text-align: center;
            padding: 20px;
            border: 1px solid #ccc;
            background-color: #fbfaf8;
        }}
        .metric-value {{
            font-size: 32px;
            font-weight: bold;
            color: #8b4513;
            display: block;
        }}
        .metric-label {{
            font-size: 14px;
            color: #666;
            margin-top: 5px;
        }}
        .charts {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin: 40px 0;
        }}
        .chart-container {{
            text-align: center;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-family: "Courier New", monospace;
            font-size: 12px;
        }}
        th, td {{
            border: 1px solid #999;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #f5f4f2;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #fbfaf8;
        }}
        .note {{
            margin: 20px 0;
            padding: 15px;
            background-color: #f9f9f9;
            border-left: 4px solid #8b4513;
            font-style: italic;
            color: #666;
        }}
    </style>
</head>
<body>
    <h1>MSI Quantification Report</h1>

    <div class="header-info">
        Microsatellite Instability Quantification Analysis<br>
        Generated {timestamp}
    </div>

    <div class="metrics">
        <div class="metric">
            <span class="metric-value">{overview['msi_rate']:.1f}%</span>
            <div class="metric-label">MSI Rate</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['msi_regions']:,}</span>
            <div class="metric-label">MSI Unstable Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['non_msi_regions_with_indels']:,}</span>
            <div class="metric-label">Uncertain Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{overview['stable_regions']:,}</span>
            <div class="metric-label">MSI Stable Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{scope['total_ms_regions']:,}</span>
            <div class="metric-label">Total MS Regions</div>
        </div>
        <div class="metric">
            <span class="metric-value">{indels['insertion_count'] + indels['deletion_count']:,}</span>
            <div class="metric-label">Perfect MSI Indels</div>
        </div>
        <div class="metric">
            <span class="metric-value">{af['mean']:.2f}</span>
            <div class="metric-label">Mean AF_max (non-zero, non-N/A)</div>
        </div>
    </div>

    <div class="note">
        <strong>Analysis Scope:</strong> {scope['analysis_note']} 
        MSI stability refers to microsatellite-specific analysis only. Other variant types (SNVs, structural variants) were not analyzed.
    </div>

    <div class="charts">
        <div class="chart-container">
            <div id="overview-chart"></div>
        </div>
        <div class="chart-container">
            <div id="simple-indel-chart"></div>
        </div>
        <div class="chart-container">
            <div id="indel-histogram"></div>
        </div>
        <div class="chart-container">
            <div id="af-chart"></div>
        </div>
        <div class="chart-container">
            <div id="motif-chart"></div>
        </div>
    </div>




    <h2>Allele Frequency Analysis</h2>
    <table>
        <thead>
            <tr>
                <th>AF_max Category</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>Variants with AF_max > 0</td>
                <td>{af['count']:,}</td>
                <td>{(af['count']/scope['perfect_msi_variants_analyzed']*100):.1f}%</td>
            </tr>
            <tr>
                <td>Variants with AF_max = 0</td>
                <td>{af['count_with_zero_af']:,}</td>
                <td>{(af['count_with_zero_af']/scope['perfect_msi_variants_analyzed']*100):.1f}%</td>
            </tr>
            <tr>
                <td>Variants with AF_max = N/A</td>
                <td>{af['count_with_na_af']:,}</td>
                <td>{(af['count_with_na_af']/scope['perfect_msi_variants_analyzed']*100):.1f}%</td>
            </tr>
            <tr>
                <td><strong>Total Perfect MSI Variants</strong></td>
                <td><strong>{scope['perfect_msi_variants_analyzed']:,}</strong></td>
                <td><strong>100.0%</strong></td>
            </tr>
        </tbody>
    </table>
    <div class="note">
        {af['mean_note']}
    </div>

    <h2>Motif Type Breakdown</h2>
    <table>
        <thead>
            <tr>
                <th>Motif Type</th>
                <th>Total Regions</th>
                <th>MSI Unstable</th>
                <th>MSI Rate</th>
                <th>Uncertain (N/A indels)</th>
                <th>MSI Stable</th>
            </tr>
        </thead>
        <tbody>"""

    for motif_type, data in sorted(msi_data["motif_breakdown"].items()):
        msi_rate = (data["msi"] / data["total"] * 100) if data["total"] > 0 else 0
        html_content += f"""
            <tr>
                <td>{motif_type}</td>
                <td>{data['total']:,}</td>
                <td>{data['msi']:,}</td>
                <td>{msi_rate:.1f}%</td>
                <td>{data['non_msi_with_indels']:,}</td>
                <td>{data['stable']:,}</td>
            </tr>"""

    html_content += f"""
        </tbody>
    </table>
    <div class="note">
        MSI Stable = no MSI indels detected in region. Uncertain = regions with N/A indels that may represent complex MSI events.
    </div>

    <script>
        vegaEmbed('#overview-chart', {overview_chart.to_json()}, {{actions: false}});
        vegaEmbed('#simple-indel-chart', {simple_indel_chart.to_json()}, {{actions: false}});
        vegaEmbed('#indel-histogram', {indel_histogram.to_json()}, {{actions: false}});
        vegaEmbed('#af-chart', {af_distribution.to_json()}, {{actions: false}});
        vegaEmbed('#motif-chart', {motif_breakdown_chart.to_json()}, {{actions: false}});
    </script>
</body>
</html>
    """

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"[MSI-ANALYSIS] HTML report saved to: {output_path}")


def create_overview_chart(msi_data):
    """Create overview pie chart of region classification"""
    overview = msi_data["overview"]

    data = [
        {
            "Category": "MSI Unstable",
            "Count": overview["msi_regions"],
            "Color": "#8b4513",
        },
        {
            "Category": "Uncertain",
            "Count": overview["non_msi_regions_with_indels"],
            "Color": "#daa520",
        },
        {
            "Category": "MSI Stable",
            "Count": overview["stable_regions"],
            "Color": "#2f4f4f",
        },
    ]

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_arc()
        .encode(
            theta=alt.Theta(field="Count", type="quantitative"),
            color=alt.Color(field="Color", type="nominal", scale=None),
            tooltip=[
                alt.Tooltip("Category:N", title="Category"),
                alt.Tooltip("Count:Q", title="Count", format=","),
            ],
        )
        .properties(width=200, height=200, title="MSI Region Classification")
    )

    return chart


def create_simple_indel_chart(msi_data):
    """Create simple histogram of insertion vs deletion counts with percentages"""
    indels = msi_data["indels"]
    total = indels["insertion_count"] + indels["deletion_count"]

    data = [
        {
            "Type": "Insertions",
            "Count": indels["insertion_count"],
            "Percentage": (indels["insertion_count"] / total * 100) if total > 0 else 0,
            "Label": (
                f"{indels['insertion_count']:,} ({(indels['insertion_count']/total*100):.1f}%)"
                if total > 0
                else "0"
            ),
            "Color": "#8b4513",
        },
        {
            "Type": "Deletions",
            "Count": indels["deletion_count"],
            "Percentage": (indels["deletion_count"] / total * 100) if total > 0 else 0,
            "Label": (
                f"{indels['deletion_count']:,} ({(indels['deletion_count']/total*100):.1f}%)"
                if total > 0
                else "0"
            ),
            "Color": "#2f4f4f",
        },
    ]

    bars = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(width=60)
        .encode(
            x=alt.X("Type:N", title="Indel Type"),
            y=alt.Y("Count:Q", title="Count"),
            color=alt.Color(field="Color", type="nominal", scale=None),
            tooltip=[
                alt.Tooltip("Type:N", title="Type"),
                alt.Tooltip("Count:Q", title="Count", format=","),
                alt.Tooltip("Percentage:Q", title="Percentage", format=".1f"),
            ],
        )
    )

    text = (
        alt.Chart(pl.DataFrame(data))
        .mark_text(align="center", baseline="bottom", dy=-5, fontSize=9)
        .encode(x=alt.X("Type:N"), y=alt.Y("Count:Q"), text=alt.Text("Label:N"))
    )

    chart = (bars + text).properties(
        width=200, height=150, title="Perfect MSI Indel Counts"
    )

    return chart


def create_indel_histogram(msi_data):
    """Create histogram of indel sizes with better hover and bins"""
    indels = msi_data["indels"]

    data = []
    for size in indels["insertion_sizes"]:
        data.append({"Size": size, "Type": "Insertion"})
    for size in indels["deletion_sizes"]:
        data.append({"Size": size, "Type": "Deletion"})

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(opacity=0.7, stroke="white", strokeWidth=0.5)
        .encode(
            x=alt.X("Size:Q", bin=alt.Bin(maxbins=15), title="Indel Size (bp)"),
            y=alt.Y("count():Q", title="Count"),
            color=alt.Color("Type:N", scale=alt.Scale(range=["#8b4513", "#2f4f4f"])),
            tooltip=[
                alt.Tooltip("Size:Q", title="Size Range (bp)", bin=alt.Bin(maxbins=15)),
                alt.Tooltip("count():Q", title="Count"),
                alt.Tooltip("Type:N", title="Type"),
            ],
        )
        .properties(width=200, height=150, title="MSI Indel Size Distribution")
    )

    return chart


def create_af_distribution_chart(msi_data):
    """Create AF_max distribution histogram (excluding 0 and N/A)"""
    af_values = msi_data["allele_frequencies"]["values"]

    if not af_values:
        return (
            alt.Chart(pl.DataFrame([{"x": 0, "y": 0}]))
            .mark_text(text="No AF data", size=14)
            .encode(x="x:Q", y="y:Q")
            .properties(width=200, height=150, title="AF_max Distribution")
        )

    data = [{"AF": af} for af in af_values]

    chart = (
        alt.Chart(pl.DataFrame(data))
        .mark_bar(color="#556b2f")
        .encode(
            x=alt.X("AF:Q", bin=alt.Bin(maxbins=15), title="AF_max"),
            y=alt.Y("count():Q", title="Count"),
            tooltip=[
                alt.Tooltip("AF:Q", title="AF_max Range", format=".2f"),
                alt.Tooltip("count():Q", title="Count"),
            ],
        )
        .properties(width=200, height=150, title="AF_max Distribution (>0 only)")
    )

    return chart


def create_motif_breakdown_chart(msi_data):
    """Create motif type MSI rate chart"""
    motif_data = []
    for motif_type, data in msi_data["motif_breakdown"].items():
        msi_rate = (data["msi"] / data["total"] * 100) if data["total"] > 0 else 0
        motif_data.append(
            {"Motif": motif_type, "MSI_Rate": msi_rate, "Total": data["total"]}
        )

    chart = (
        alt.Chart(pl.DataFrame(motif_data))
        .mark_bar(color="#8b4513")
        .encode(
            x=alt.X("Motif:N", title="Motif Type", sort="-y"),
            y=alt.Y(
                "MSI_Rate:Q", title="MSI Rate (%)", scale=alt.Scale(domain=[0, 40])
            ),
            tooltip=[
                alt.Tooltip("Motif:N", title="Motif Type"),
                alt.Tooltip("MSI_Rate:Q", title="MSI Rate (%)", format=".1f"),
                alt.Tooltip("Total:Q", title="Total Regions", format=","),
            ],
        )
        .properties(width=200, height=150, title="MSI Rate by Motif Type")
    )

    return chart


################################################################################################


def main():
    """Main function to execute the MSI analysis."""

    args = parse_args()

    initialize_debug_log(args.debug_log, args)

    print("Loading VCF variants...")
    variants, filter_stats = load_vcf_variants(args.vcf)

    print("Loading BED regions...")
    regions = load_bed_regions(args.bed)

    print(f"Running intersection analysis with {args.threads} threads...")
    results, total_regions_loaded, unprocessed_count, merged_count = (
        bin_based_threaded_intersection(variants, regions, args.threads)
    )

    output_data = {
        "analysis_info": {
            "timestamp": datetime.now().isoformat(),
            "vcf_file": args.vcf,
            "bed_file": args.bed,
            "threads_used": args.threads,
        },
        "region_results": results,
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=2, default=str)

    print(f"MSI analysis complete! Results saved to {args.output}")

    write_complete_debug_log(
        args.debug_log,
        variants,
        regions,
        results,
        total_regions_loaded,
        args,
        filter_stats,
        unprocessed_count,
        merged_count,
    )

    msi_data = analyze_msi_results_post_quantification(results)
    generate_msi_html_report(msi_data, args.html_report)

    print(f"Debug log: {args.debug_log}")
    print(f"HTML report: {args.html_report}")



if __name__ == "__main__":
    main()
