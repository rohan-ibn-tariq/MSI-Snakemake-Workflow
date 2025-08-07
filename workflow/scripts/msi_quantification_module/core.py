"""
MSI Quantification Module - Core Analysis Functions
Core MSI analysis logic, data loading, and intersection processing.
"""

import multiprocessing
from collections import defaultdict

import pysam
from binning import assign_bin
from scipy.stats import chisquare
from scipy.stats import false_discovery_control


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
                reference_length = repeat_count * motif_length
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
                    "reference_length": reference_length,
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
        "prob_present": variant_data["prob_present"],
        "prob_absent": variant_data["prob_absent"],
        "prob_artifact": variant_data["prob_artifact"],
        "sample_afs": variant_data["sample_afs"],
    }


# def analyze_variant_in_region(variant_data, region_data):
#     """Analyze a variant within a microsatellite region - DEBUG VERSION"""
    
#     # DEBUG: Track the problematic region
#     debug_region = "10518525" in region_data.get("region_id", "")
#     if debug_region:
#         print(f"    ANALYZING: {variant_data['pos']} {variant_data['ref']}→{variant_data['alt']}")

#     if not check_variant_overlaps_region(variant_data, region_data):
#         if debug_region:
#             print(f"      REJECTED: No overlap")
#         return None

#     ref_seq = variant_data["ref"]
#     alt_seq = variant_data["alt"]

#     if debug_region:
#         print(f"      ref_seq='{ref_seq}', alt_seq='{alt_seq}'")

#     anchor_len = 0
#     min_len = min(len(ref_seq), len(alt_seq))

#     for i in range(min_len):
#         if ref_seq[i].upper() == alt_seq[i].upper():
#             anchor_len += 1
#         else:
#             break

#     if debug_region:
#         print(f"      anchor_len={anchor_len}, min_len={min_len}")

#     if anchor_len == 0:
#         if debug_region:
#             print(f"      REJECTED: No anchor base found")
#         print(f"[WARNING] No anchor base found for variant {variant_data['variant_id']}")
#         return None

#     if variant_data["svlen"] > 0:
#         changed_seq = alt_seq[anchor_len:]
#     else:
#         changed_seq = ref_seq[anchor_len:]

#     if debug_region:
#         print(f"      svlen={variant_data['svlen']}, changed_seq='{changed_seq}'")

#     if not changed_seq:
#         if debug_region:
#             print(f"      REJECTED: No changed sequence")
#         print(f"[WARNING] No changed sequence for variant {variant_data['variant_id']}")
#         return None

#     motif = region_data["motif"]
#     is_perfect = is_perfect_repeat(changed_seq, motif)
#     repeat_status = "perfect" if is_perfect else "N/A"

#     if debug_region:
#         print(f"      motif='{motif}', is_perfect={is_perfect}, repeat_status='{repeat_status}'")
#         print(f"      ACCEPTED: Creating analysis result")

#     return {
#         "variant_id": variant_data["variant_id"],
#         "vcf_id": variant_data.get("vcf_id"),
#         "region_id": region_data["region_id"],
#         "motif": motif,
#         "motif_type": region_data["motif_type"],
#         "variant_type": variant_data["variant_type"],
#         "svlen": variant_data["svlen"],
#         "changed_sequence": changed_seq,
#         "repeat_status": repeat_status,
#         "af_mean": variant_data["af_mean"],
#         "af_max": variant_data["af_max"],
#     }


def calculate_chi_squared_for_perfect_variants(region, perfect_variants):
    """
    Calculate chi-squared test for MSI instability at a single region
    """
    ref_length = region["reference_length"]

    observed = {}
    for variant in perfect_variants:
        alt_length = ref_length + variant["svlen"]
        observed[alt_length] = observed.get(alt_length, 0) + 1

    if ref_length not in observed:
        observed[ref_length] = 0

    all_lengths = sorted(observed.keys())
    total = sum(observed.values())

    observed_counts = [observed[l] + 0.5 for l in all_lengths]

    expected_counts = []
    for l in all_lengths:
        if l == ref_length:
            expected_counts.append(total + 0.5)
        else:
            expected_counts.append(0.5)

    chi2_stat, p_value = chisquare(f_obs=observed_counts, f_exp=expected_counts)

    return chi2_stat, p_value


# def add_chi_squared_classification(results):
#     """
#     Add chi-squared MSI classification to existing results
#     """
#     unstable_regions = 0
#     testable_regions = 0
#     total_regions = len(results)

#     for region in results:
#         if region.get("num_perfect_repeats", 0) > 0:
#             testable_regions += 1

#             perfect_variants = [v for v in region["variants"]
#                               if v["repeat_status"] == "perfect"]

#             chi2_stat, p_value = calculate_chi_squared_for_perfect_variants(region, perfect_variants)

#             ###DEBUG
#             # ADD DEBUG PRINT
#             if testable_regions <= 5:  # Print first 5 for debugging
#                 print(f"Region {testable_regions}: chi2={chi2_stat:.3f}, p={p_value:.6f}, variants={len(perfect_variants)}")
#             ###END

#             region["is_unstable"] = p_value < 0.05
#             region["chi2_stat"] = chi2_stat
#             region["p_value"] = p_value

#             if region["is_unstable"]:
#                 unstable_regions += 1
#         else:
#             region["is_unstable"] = False
#             region["chi2_stat"] = None
#             region["p_value"] = None

#     ###DEBUG

#     # Count regions by p-value ranges
#     p_ranges = {"< 0.001": 0, "0.001-0.01": 0, "0.01-0.05": 0, "> 0.05": 0}
#     for region in results:
#         if region.get("p_value") is not None:
#             p = region["p_value"]
#             if p < 0.001:
#                 p_ranges["< 0.001"] += 1
#             elif p < 0.01:
#                 p_ranges["0.001-0.01"] += 1
#             elif p < 0.05:
#                 p_ranges["0.01-0.05"] += 1
#             else:
#                 p_ranges["> 0.05"] += 1

#     print("P-value distribution:", p_ranges)
#     ###END
#     msi_score = (unstable_regions / total_regions) * 100 if total_regions > 0 else 0

#     if msi_score >= 3.5:
#         msi_status = "MSI-High"
#     else:
#         msi_status = "MSS"

#     return {
#         "msi_score": msi_score,
#         "msi_status": msi_status,
#         "unstable_regions": unstable_regions,
#         "testable_regions": testable_regions,
#         "total_regions": total_regions,
#         "threshold_used": 3.5
#     }


def add_chi_squared_classification(results):
    """
    Add chi-squared MSI classification with FDR correction
    """
    testable_regions = 0
    total_regions = len(results)

    p_values = []
    regions_with_tests = []

    for region in results:
        if region.get("num_perfect_repeats", 0) > 0:
            testable_regions += 1

            perfect_variants = [
                v for v in region["variants"] if v["repeat_status"] == "perfect"
            ]

            chi2_stat, p_value = calculate_chi_squared_for_perfect_variants(
                region, perfect_variants
            )

            region["chi2_stat"] = chi2_stat
            region["p_value"] = p_value

            p_values.append(p_value)
            regions_with_tests.append(region)
        else:
            region["is_unstable"] = False
            region["chi2_stat"] = None
            region["p_value"] = None

    ###DEBUG START
    print(f"Total testable regions: {testable_regions}")
    print(f"P-values collected: {len(p_values)}")
    if p_values:
        print(f"Original p-values (first 10): {p_values[:10]}")
        print(f"Min original p-value: {min(p_values)}")
        print(f"Max original p-value: {max(p_values)}")
        print(f"Unique p-values: {len(set(p_values))}")
    ###DEBUG END

    # Apply FDR correction (Benjamini-Hochberg)
    unstable_regions = 0
    if p_values:
        corrected_p_values = false_discovery_control(p_values, method="bh")

        ###DEBUG START
        print(f"Corrected p-values (first 10): {corrected_p_values[:10]}")
        print(f"Min corrected p-value: {min(corrected_p_values)}")
        print(f"Max corrected p-value: {max(corrected_p_values)}")

        significant_001 = sum(1 for p in corrected_p_values if p < 0.01)
        significant_005 = sum(1 for p in corrected_p_values if p < 0.05)
        significant_010 = sum(1 for p in corrected_p_values if p < 0.10)

        print(f"Corrected p-values < 0.01: {significant_001}")
        print(f"Corrected p-values < 0.05: {significant_005}")
        print(f"Corrected p-values < 0.10: {significant_010}")
        ###DEBUG END

        for i, region in enumerate(regions_with_tests):
            region["p_value_corrected"] = corrected_p_values[i]
            region["is_unstable"] = corrected_p_values[i] < 0.05

            if region["is_unstable"]:
                unstable_regions += 1

    ###DEBUG START
    print(f"Final unstable regions: {unstable_regions}")
    ###DEBUG END

    msi_score = (unstable_regions / total_regions) * 100 if total_regions > 0 else 0

    if msi_score >= 3.5:
        msi_status = "MSI-High"
    else:
        msi_status = "MSS"

    return {
        "msi_score": msi_score,
        "msi_status": msi_status,
        "unstable_regions": unstable_regions,
        "testable_regions": testable_regions,
        "total_regions": total_regions,
        "threshold_used": 3.5,
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
            "reference_length": region_data["reference_length"],
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

    msi_classification = add_chi_squared_classification(results)
    print(msi_classification)

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
        "msi_classification": msi_classification,
    }
