"""
MSI Quantification Module - Dynamic Programming Analysis

Implements probabilistic MSI analysis using dynamic programming approach.
Calculates MSI probability distributions and uncertainty quantification across AF thresholds.

Outputs:
1. AF Evolution Analysis: MSI scores across allele frequency thresholds per sample
2. MSI Probability Distributions: Complete probability tables for AF 0.0
3. MAP-based Uncertainty: Statistical confidence intervals around MSI estimates

Features:
- Trusts upstream AF validation (0.0-1.0 range)
- Missing data handling (-1 for null AF values)
- MAP-centered uncertainty calculation
- Memory-efficient lazy evaluation for large datasets

TODO: Consider upstream DP integration for single-pass processing
"""

from typing import Dict, List, Optional, Union

import pysam


#####################################################################################
################## DATA QUALITY VALIDATION & IMPUTATION FUNCTIONS ###################
#####################################################################################


def validate_af_value(af_value: Union[float, int, str, None]) -> Union[float, int]:
    """
    Normalize AF value with upstream validation trust.

    Args:
        af_value: AF value from VCF (any type or None)

    Returns:
        float (0.0-1.0) for valid values, -1 for missing/null
    """
    if af_value is None:
        return -1

    return float(af_value)


def process_variant_probabilities_and_af(variant: Dict) -> None:
    """
    Process variant probabilities and AF values for DP analysis.

    Side Effects:
        - Adds variant["dp_data"] field with DP-ready probabilities and normalized AF values
    """
    pp = variant.get("prob_present")
    pa = variant.get("prob_absent")
    art = variant.get("prob_artifact")

    af_by_sample = {}
    for sample_name, af_value in variant.get("sample_afs", {}).items():
        af_by_sample[sample_name] = validate_af_value(af_value)

    variant["dp_data"] = {
        "p_present": pp,
        "p_absent": pa + art,
        "af_by_sample": af_by_sample,
    }


def extract_af_data_for_variant(
    variant: Dict, sample_list: List[str]
) -> Dict[str, Union[float, int]]:
    """
    Extract AF values for N/A variants in uncertain region analysis.

    Args:
        variant: N/A variant dictionary containing sample_afs field
        sample_list: Complete list of sample names from VCF header

    Returns:
        Dict mapping sample names to AF values:
            - float (0.0-1.0): Valid AF value
            - -1: Missing/null AF or sample not present in variant

    TODO: Consider upstream consolidation - move all AF processing to core analysis
          to eliminate preparation phase and reduce redundant AF validation calls.
    """
    sample_afs = variant.get("sample_afs", {})
    return {name: validate_af_value(sample_afs.get(name)) for name in sample_list}


#####################################################################################
################## DATA PREPARATION #################################################
#####################################################################################


def prepare_variants_for_dp(results, vcf_file_path, imputation_method="uniform"):
    """
    Data Preparation Function for Dynamic Programming Analysis.

    #TODO: ARCHITECTURE OPTIMIZATION - Move data preparation upstream
    # Current: Core analysis → DP preparation → DP analysis (3 stages)
    # Alternative: Integrate DP data preparation into analyze_variant_in_region()
    # or bin based intersection upstream in core analysis.
    # Benefits: Single-pass processing, eliminate separate preparation stage
    # Trade-off: Core analysis becomes heavier vs current clean separation
    # Note: Requires passing sample_list to core analysis
    # State: Current separation provides good modularity, optimize only if bottleneck confirmed
    """

    print("[DP-ENTRY] Starting data preparation for dynamic programming analysis")
    print(f"[DP-ENTRY] Processing {len(results)} regions")

    vcf = pysam.VariantFile(vcf_file_path)
    sample_list = list(vcf.header.samples)
    vcf.close()
    print(f"[DP-PREP] Found samples in VCF: {sample_list}")

    af_thresholds = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1]

    total_uncertain_regions = 0

    af_data_by_sample = {}
    for sample_name in sample_list:
        af_data_by_sample[sample_name] = {}
        for af_threshold in af_thresholds:
            af_data_by_sample[sample_name][af_threshold] = {
                "uncertain_regions": 0,
                "variant_count": 0,
            }

    # Extract variants and process in single pass
    region_variants = {}
    total_variants = 0

    for region in results:
        region_id = region.get("region_id")
        variants = region.get("variants", [])

        perfect_variants = [v for v in variants if v.get("repeat_status") == "perfect"]

        na_variants = [v for v in variants if v.get("repeat_status") == "N/A"]
        if na_variants and not perfect_variants:
            total_uncertain_regions += 1

            # Per-sample counting for AF evolution of uncertain regions
            for sample_name in sample_list:
                for af_threshold in af_thresholds:
                    region_has_variants_at_af_for_sample = False
                    for variant in na_variants:
                        af_data = extract_af_data_for_variant(variant, sample_list)
                        # TODO: IN REFACTORING FIX
                        sample_af = af_data[sample_name]
                        if (af_threshold == -1 and sample_af == -1) or (
                            af_threshold >= 0 and sample_af >= af_threshold
                        ):
                            region_has_variants_at_af_for_sample = True
                            break

                    if region_has_variants_at_af_for_sample:
                        af_data_by_sample[sample_name][af_threshold][
                            "uncertain_regions"
                        ] += 1

        if perfect_variants:
            # Process each variant in this region
            for variant in perfect_variants:
                # Calculate DP probabilities and AF values
                process_variant_probabilities_and_af(variant)

                # Collect sample names
                af_by_sample = variant["dp_data"]["af_by_sample"]
                for sample_name in sample_list:
                    if sample_name not in af_by_sample:
                        af_by_sample[sample_name] = -1

                    sample_af = af_by_sample[sample_name]

                    for af_threshold in af_thresholds:
                        if (af_threshold == -1 and sample_af == -1) or (
                            af_threshold >= 0 and sample_af >= af_threshold
                        ):
                            af_data_by_sample[sample_name][af_threshold][
                                "variant_count"
                            ] += 1

                total_variants += 1

            # Store processed variants
            region_variants[region_id] = perfect_variants

    # Handle edge case
    if not region_variants:
        return {
            "error": "No variants found",
            "step": "data_preparation_failed",
            "total_regions": len(results),
        }

    # Create final structure
    sample_list = sorted(sample_list)

    dp_ready = {
        "af_thresholds": af_thresholds,
        "samples": sample_list,
        "regions": region_variants,
        "total_variants": total_variants,
        "total_regions": len(region_variants),
        "ready_for_dp": True,
        "af_data_by_sample": af_data_by_sample,
        "total_uncertain_regions": total_uncertain_regions,
    }

    print(
        f"[DP-PREP] Processed {total_variants} variants in {len(region_variants)} regions"
    )
    print(f"[DP-PREP] Found {total_uncertain_regions} uncertain regions")
    print(f"[DP-PREP] AF thresholds: {af_thresholds}")
    print("[DP-PREP] Data Preparation for DP complete")

    return dp_ready


#####################################################################################
################## DYNAMIC PROGRAMMING CORE  ########################################
#####################################################################################


def run_msi_dp(variants_with_probabilities):
    """
    Execute dynamic programming algorithm for MSI variant probability distribution.

    Implements matrix-based approach where:
    - Matrix has n columns (variants 0 to n-1) and n+1 rows (MSI counts 0 to n)
    - Our implementation stores row vectors of length n+1
    - Rows represent possible MSI counts (0 to n)
    - Each cell contains P(exactly i MSI variants using first k variants)

    Algorithm:
    - Initialize column 0: m[0,0] = p_absent_0, m[1,0] = p_present_0
    - Recurrence: m[i,k] = m[i,k-1] × p_absent_k + m[i-1,k-1] × p_present_k

    Args:
        variants_with_probabilities (List[Dict]): Variants with dp_data containing
            p_present and p_absent probabilities from 4-step imputation

    Returns:
        List[float]: Probability distribution [P(0 MSI), P(1 MSI), ..., P(n MSI)]
        where sum equals 1.0 and each P(i) represents probability of exactly
        i variants being MSI in the region

    Note:
        Uses space-optimized approach storing only previous column instead of full matrix.
        p_present + p_absent should equal 1.0 from imputation step.
        - PHRED probability conversion has ~0.005 tolerance for floating-point precision
    """
    n = len(variants_with_probabilities)

    # Base case: no variants means P(0 MSI) = 100%
    if n == 0:
        return [1.0]

    # Initialize probability vector: [P(0), P(1), P(2), ..., P(n)]
    prev_col = [0.0] * (n + 1)

    # Column 0: Initialize with first variant probabilities
    p_absent_0 = variants_with_probabilities[0]["dp_data"]["p_absent"]
    p_present_0 = variants_with_probabilities[0]["dp_data"]["p_present"]

    prev_col[0] = p_absent_0  # P(0 MSI) = first variant absent
    prev_col[1] = p_present_0  # P(1 MSI) = first variant present

    # Process remaining variants using recurrence relation
    for k in range(1, n):
        p_absent_k = variants_with_probabilities[k]["dp_data"]["p_absent"]
        p_present_k = variants_with_probabilities[k]["dp_data"]["p_present"]

        curr_col = [0.0] * (n + 1)

        # Base case: P(0 MSI) = previous P(0 MSI) × current variant absent
        curr_col[0] = prev_col[0] * p_absent_k

        # Recurrence: P(exactly i MSI) has two paths:
        # Path 1: Had i MSI, current variant absent → still i MSI
        # Path 2: Had i-1 MSI, current variant present → now i MSI
        for i in range(1, k + 2):  # Up to k+1 total MSI variants possible
            curr_col[i] = prev_col[i] * p_absent_k + prev_col[i - 1] * p_present_k

        prev_col = curr_col  # Update for next iteration

    return prev_col  # Final distribution: [P(0), P(1), ..., P(n)]


#####################################################################################
############################## HELPERS FOR DP ANALYSIS ##############################
#####################################################################################
def classify_msi_status(msi_score: float, msi_high_threshold: float = 3.5) -> str:
    """
    Classify MSI status using MSIsensor standard threshold.

    Args:
        msi_score (float): MSI score as percentage (0.0-100.0)
        msi_high_threshold (float, optional): Threshold for MSI-High classification.
            Defaults to 3.5 (MSIsensor standard).

    Returns:
        str: MSI classification ("MSI-High" or "MSS")

    Raises:
        ValueError: If msi_score is negative or msi_high_threshold is non-positive
    """
    if msi_score < 0:
        raise ValueError(f"MSI score must be non-negative, got {msi_score}")
    if msi_high_threshold <= 0:
        raise ValueError(f"MSI threshold must be positive, got {msi_high_threshold}")

    return "MSI-High" if msi_score >= msi_high_threshold else "MSS"


#####################################################################################
################## DP ANALYSIS FUNCTIONS ############################################
#####################################################################################


def calculate_msi_metrics_for_regions(
    regions_dict: Dict[str, List],
    total_regions: int,
    uncertain_regions: int,
    msi_high_threshold: float = 3.5,
    af_threshold: Optional[float] = None,
) -> Dict:
    """
    Calculate MSI metrics based on DP distributions per region.
    Returns MSI data per sample (prob. distribution, MSI score, MAP score),
    plus summary counts for reporting.
    """

    # Compute per-region instability probabilities
    region_probabilities = []
    for _, variants in regions_dict.items():
        p_zero = 1.0
        for variant in variants:
            p_zero *= variant["dp_data"]["p_absent"]
        region_probabilities.append(1 - p_zero)

    regions_with_variants = len(region_probabilities)

    # Build distribution of unstable counts using DP
    if region_probabilities:
        region_variants_for_dp = [
            {"dp_data": {"p_present": p, "p_absent": 1 - p}}
            for p in region_probabilities
        ]
        distribution = run_msi_dp(region_variants_for_dp)
    else:
        distribution = [1.0]  # no regions → 0 unstable with prob 1

    # Prepare MSI data per sample count for AF >= 0.0 only
    msi_data = None
    if af_threshold == 0.0:
        msi_data = {}
        for k, prob in enumerate(distribution):
            msi_score = (k / total_regions) * 100 if total_regions > 0 else 0.0
            msi_data[k] = {
                "probability": prob,
                "msi_score": round(msi_score, 2),
            }

    # MAP estimate
    k_map = max(range(len(distribution)), key=lambda i: distribution[i])
    msi_score_map = (k_map / total_regions) * 100 if total_regions > 0 else 0.0

    # MAP-based uncertainty
    map_variance = sum((k - k_map) ** 2 * prob for k, prob in enumerate(distribution))
    map_std = map_variance**0.5

    # Convert to MSI score uncertainty range
    uncertainty_msi_lower = max(0, (k_map - map_std) / total_regions * 100)
    uncertainty_msi_upper = min(100, (k_map + map_std) / total_regions * 100)
    uncertainty_range = [uncertainty_msi_lower, uncertainty_msi_upper]

    base_result = {
        # Counts / reporting
        "total_regions": total_regions,
        "regions_with_variants": regions_with_variants,
        "uncertain_regions": uncertain_regions,
        # MSI per-sample data
        "k_map": k_map,
        "msi_score_map": round(msi_score_map, 2),
        "msi_status_map": classify_msi_status(msi_score_map, msi_high_threshold),
        "uncertainty_range": uncertainty_range,
        "map_std_dev": round(map_std, 3),
        # Analysis parameters
        "analysis_parameters": {
            "msi_high_threshold": msi_high_threshold,
        },
    }

    if msi_data is not None:
        base_result["msi_data"] = msi_data

    return base_result


def filter_variants_by_af_and_sample(
    variants: List[Dict], af_threshold: float, sample_name: str
) -> List[Dict]:
    """
    Filter variants by AF threshold for a specific sample.

    Args:
        variants: List of variants in a region
        af_threshold: Minimum AF threshold
        sample_name: Sample to filter for

    Returns:
        List of variants that meet AF threshold for this sample

    #TODO: AF PRE-COMPUTATION OPTIMIZATION
    # Current: This function called repeatedly in nested loops (samples × thresholds × regions)
    # Results in ~n+ filtering operations for same data
    # Alternative: Pre-compute all AF-filtered combinations in prepare_variants_for_dp()
    # Trade-off: ~x% performance gain vs z×-y× memory increase
    # State: Current lazy evaluation chosen for memory efficiency and scalability
    """
    filtered_variants = []
    additional_uncertain = 0

    variants_meeting_threshold = []
    variants_with_missing_af = []

    for variant in variants:
        sample_af = variant["dp_data"]["af_by_sample"][sample_name]
        if sample_af >= af_threshold:
            filtered_variants.append(variant)
            variants_meeting_threshold.append(variant)
        elif sample_af == -1:
            variants_with_missing_af.append(variant)

    # If region has missing AF variants but no variants meeting threshold
    if not variants_meeting_threshold and variants_with_missing_af:
        additional_uncertain = 1

    return filtered_variants, additional_uncertain


def run_af_evolution_analysis(
    dp_ready_data: Dict,
    total_ms_regions: int,
    msi_high_threshold: float = 3.5,
) -> Dict:
    """
    Run AF-based MSI evolution analysis across allele frequency thresholds.

    Analyzes MSI scores at different AF thresholds for each sample, providing
    insights into tumor evolution from early (high AF) to late (low AF) mutations.
    Uses the same DP analysis as regional analysis but applied to AF-filtered variants.

    Args:
        dp_ready_data (Dict): Output from prepare_variants_for_dp()
        total_ms_regions (int): Total MS regions from BED file
        unstable_threshold (float): P(≥1 MSI) threshold for calling region unstable (default: 0.5)
        msi_high_threshold (float): MSI score threshold for MSI-High classification (default: 3.5)

    Returns:
        Dict: AF evolution results grouped by sample and AF threshold
    """
    print("[AF-EVOLUTION] Starting simplified AF evolution analysis")

    results = {}
    samples = dp_ready_data["samples"]
    af_thresholds = dp_ready_data["af_thresholds"]

    valid_af_thresholds = [t for t in af_thresholds if t >= 0.0]

    for sample_name in samples:
        print(f"[AF-EVOLUTION] Analyzing sample: {sample_name}")
        sample_results = {}

        for af_threshold in valid_af_thresholds:
            filtered_regions = {}
            total_additional_uncertain = 0

            for region_id, variants in dp_ready_data["regions"].items():
                filtered_variants, uncertain_adjustment = (
                    filter_variants_by_af_and_sample(
                        variants, af_threshold, sample_name
                    )
                )
                total_additional_uncertain += uncertain_adjustment

                if filtered_variants:
                    filtered_regions[region_id] = filtered_variants

            af_data = dp_ready_data["af_data_by_sample"][sample_name][af_threshold]
            uncertain_regions = (
                af_data["uncertain_regions"] + total_additional_uncertain
            )
            # total_variants = af_data["variant_count"]

            metrics = calculate_msi_metrics_for_regions(
                filtered_regions,
                total_ms_regions,
                uncertain_regions,
                # total_variants,
                # sum(len(v) for v in filtered_regions.values()), #TODO: CONFIRM REQUIRED AND RIGHT OR FILTERED REQUIRED
                msi_high_threshold,
                af_threshold=af_threshold,
            )

            print(
                f"[AF-EVOLUTION] {sample_name} | AF={af_threshold} | "
                f"k_map={metrics['k_map']} | "
                f"MSI%={metrics['msi_score_map']} | "
                f"Regions={metrics['regions_with_variants']}/{metrics['total_regions']} "
                f"(uncertain={uncertain_regions})"
            )

            sample_results[f"af_{af_threshold}"] = metrics

        results[sample_name] = sample_results

    return results
