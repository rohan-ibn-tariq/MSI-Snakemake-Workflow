"""
MSI Quantification Module - Dynamic Programming Analysis

Implements probabilistic MSI analysis using dynamic programming approach.
Outputs:
1. Genome Wide Regional MSI Analysis
2. Sample Clustered Genome-Wide AF Evolution Analysis
Provides uncertainty quantification for MSI scores, unstable regions, and variants.
"""

from typing import Dict, List, Optional, Tuple, Union

import pysam

#####################################################################################
################## DATA QUALITY VALIDATION & IMPUTATION FUNCTIONS ###################
#####################################################################################

def validate_af_value(
    af_value: Union[float, int, str, None], sample_name: str
) -> Tuple[Union[float, int], List[str]]:
    """
    Validate and normalize AF value for a specific sample.

    Args:
        af_value: Allele frequency value from VCF (float, int, str, or None)
        sample_name: Sample identifier for error reporting

    Returns:
        Tuple of (normalized_af, error_list)
        - normalized_af: float (0.0-1.0), -1 (missing), or -2 (invalid)
        - error_list: List of validation errors
    """
    if af_value is None:
        return -1, []  # Missing/null AF (biological uncertainty)

    try:
        af_float = float(af_value)
        if af_float < 0:
            return -2, [f"{sample_name}_negative_af_{af_float}"]  # Invalid data
        if af_float > 1.0:
            return -2, [f"{sample_name}_af_exceeds_1.0_{af_float}"]  # Invalid data
        return af_float, []
    except (ValueError, TypeError):
        return -2, [f"{sample_name}_invalid_af_type"]  # Invalid data


# pylint: disable=too-many-locals,too-many-branches,too-many-statements
def process_variant_probabilities_and_af(variant: Dict) -> None:
    """
    Process variant probabilities and allele frequencies for DP analysis.
    Creates dp_data with consolidated probabilities and validated AF values.

    Side Effects:
        - Adds variant["dp_data"] field with results
        - Validates and normalizes all AF values with error tracking
        - Creates detailed trail for reproducibility
    
    Note: Variants with missing probabilities are filtered out by core analysis.
    """

    pp = variant.get("prob_present")
    pa = variant.get("prob_absent")
    art = variant.get("prob_artifact")

    p_present = pp
    p_absent = pa + art

    audit_trail = {
        "af_processing": {
            "missing_af_samples": [],
            "af_conversion": {},
            "af_errors": [],
        },
    }

    # Process AF values for each sample with validation and audit tracking
    sample_afs = variant.get("sample_afs", {})
    af_by_sample = {}

    for sample_name, af_value in sample_afs.items():
        normalized_af, af_errors = validate_af_value(af_value, sample_name)

        af_by_sample[sample_name] = normalized_af
        audit_trail["af_processing"]["af_errors"].extend(af_errors)

        if normalized_af == -1:
            audit_trail["af_processing"]["missing_af_samples"].append(sample_name)
            if af_value is None:
                audit_trail["af_processing"]["af_conversion"][
                    sample_name
                ] = "None -> -1"
            else:
                audit_trail["af_processing"]["af_conversion"][
                    sample_name
                ] = f"{af_value} -> -1 (error)"
        else:
            audit_trail["af_processing"]["af_conversion"][
                sample_name
            ] = f"{af_value} -> {normalized_af}"

    variant["dp_data"] = {
        "p_present": p_present,
        "p_absent": p_absent,
        "af_by_sample": af_by_sample,
        "audit_trail": audit_trail,
    }


def extract_af_data_for_variant(
    variant: Dict, sample_list: List[str]
) -> Dict[str, Union[float, int]]:
    """
    Extract and normalize AF data for N/A variants in uncertain regions.

    Used only for N/A variants during AF evolution analysis of uncertain regions.
    Provides lightweight AF extraction without the full DP processing and audit
    trails used for perfect variants in apply_4_step_imputation_to_variant.

    Args:
        variant (Dict): N/A variant dictionary containing sample_afs field
        sample_list (List[str]): Complete list of sample names from VCF header

    Returns:
        Dict[str, Union[float, int]]: AF values per sample where:
            - float (0.0-1.0): Valid AF value
            - -1: Missing/null AF (biological uncertainty)
            - -2: Invalid AF data (technical error)

    Note:
        Only called for N/A variants to count uncertain regions in AF evolution.

    #TODO: Refactor to process all variants (perfect + N/A) in single loop rather
        than separate processing. Current separation exists for time-saving
        during development but creates code duplication.
    """

    sample_afs = variant.get("sample_afs", {})
    af_by_sample = {}

    for sample_name, af_value in sample_afs.items():
        normalized_af, _ = validate_af_value(af_value, sample_name)
        af_by_sample[sample_name] = normalized_af

    for sample_name in sample_list:
        if sample_name not in af_by_sample:
            af_by_sample[sample_name] = -1

    return af_by_sample


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

    af_thresholds = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1, -2]

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
                        if (
                            (af_threshold == -1 and sample_af == -1)
                            or (af_threshold == -2 and sample_af == -2)
                            or (af_threshold >= 0 and sample_af >= af_threshold)
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
                        if (
                            (af_threshold == -1 and sample_af == -1)
                            or (af_threshold == -2 and sample_af == -2)
                            or (af_threshold >= 0 and sample_af >= af_threshold)
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
################## STATISTICAL HELPERS FOR DP ANALYSIS ##############################
#####################################################################################


def calculate_expected_value(distribution: List[float]) -> float:
    """
    Calculate expected number of MSI variants from DP probability distribution.

    Mathematical foundation: E[X] = Σ(i × P(X = i)) for i = 0 to n

    Args:
        distribution (List[float]): Probability distribution from run_msi_dp()
            where distribution[i] = P(exactly i MSI variants)

    Returns:
        float: Expected number of MSI variants

    Raises:
        ValueError: If distribution doesn't sum to 1.0 (within 0.005 tolerance)

    Note:
        Distribution from run_msi_dp() should sum to 1.0 (valid probabilities).
        Uses same 0.005 tolerance as PHRED probability conversion for consistency.

    #TODO: Fix the tolerance check to be more robust.
    """
    # Validation with PHRED tolerance
    total = sum(distribution)
    if abs(total - 1.0) > 0.005:  # Same tolerance as PHRED conversion
        raise ValueError(
            f"Invalid probability distribution: sums to {total:.6f}, expected 1.0 ± 0.005"
        )

    return sum(i * prob for i, prob in enumerate(distribution))


def calculate_expected_and_variance(distribution: List[float]) -> Tuple[float, float]:
    """
    Calculate expected value and variance from DP probability distribution.
        
    Mathematical foundation: 
    - E[X] = Σ(i × P(X = i)) for i = 0 to n
    - Var(X) = Σ((i - E[X])² × P(X = i)) for i = 0 to n
    
    Args:
        distribution: Probability distribution from run_msi_dp()
        
    Returns:
        Tuple of (expected_value, variance)

    Note:
        Used for DP uncertainty propagation in expected unstable regions and
        expected MSI variants. Distribution validation handled by calculate_expected_value.
    """
    expected = calculate_expected_value(distribution)
    variance = sum((i - expected) ** 2 * prob for i, prob in enumerate(distribution))
    return expected, variance


def calculate_p_unstable(distribution: List[float]) -> float:
    """
    Calculate probability of instability: P(≥1 MSI variant).

    Mathematical foundation: P(≥1) = 1 - P(0) = 1 - distribution[0]
    This represents the probability that at least one variant is MSI.

    Args:
        distribution (List[float]): Probability distribution from run_msi_dp()
            where distribution[0] = P(exactly 0 MSI variants)

    Returns:
        float: Probability of having ≥1 MSI variant (0.0 to 1.0)

    Note:
        Used for determining unstable regions in regional MSI analysis.
        Compared against unstable_threshold (default 0.5) for classification.
    """
    if not distribution:
        return 0.0  # No variants = no instability
    return 1.0 - distribution[0]


def propagate_uncertainties(
    variances: List[float], total_regions: int
) -> Tuple[float, float]:
    """
    Propagate statistical uncertainties for MSI score calculation using error propagation theory.

    MATHEMATICAL FOUNDATION:
    For independent random variables X₁, X₂, ..., Xₙ:
    - Variance additivity: Var(X₁ + X₂ + ... + Xₙ) = Var(X₁) + Var(X₂) + ... + Var(Xₙ)
    - Standard deviation: σ(sum) = √(σ₁² + σ₂² + ... + σₙ²)

    APPLICATION TO MSI ANALYSIS:
    - Each region contributes: variance to total unstable count variance
    - MSI score = (total unstable / total regions) × 100
    - MSI uncertainty = (sqrt(total_variance) / total regions) × 100

    Args:
        variances (List[float]): List of variances from individual regions
                            obtained from calculate_variance() applied to DP distributions
        total_regions (int): Total microsatellite regions in genome (denominator for MSI score)

    Returns:
        Tuple[float, float]: (uncertainty_in_count, msi_uncertainty_percentage)
            - uncertainty_in_count: Standard deviation of total unstable region count
            - msi_uncertainty_percentage: MSI score uncertainty in percentage points

    Raises:
        ValueError: If total_regions <= 0

    Note:
        Used for DP uncertainty propagation in regional MSI analysis to quantify
        statistical precision of MSI score estimates.
    """
    if total_regions <= 0:
        raise ValueError(f"Total_regions must be positive, got {total_regions}")

    if not variances:
        return 0.0, 0.0

    # Validate variance values
    for i, variance in enumerate(variances):
        if variance < 0:
            raise ValueError(f"Variance[{i}] must be non-negative, got {variance}")

    # Sum of variances
    total_variance = sum(variances)

    # Standard deviation of total unstable count
    uncertainty_in_count = total_variance**0.5

    # Convert to MSI score uncertainty (error propagation for division)
    # MSI uncertainty = σ_X / N
    msi_uncertainty = (uncertainty_in_count / total_regions) * 100

    return uncertainty_in_count, msi_uncertainty


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
    variants_deterministic: int,
    unstable_threshold: float = 0.5,
    msi_high_threshold: float = 3.5,
) -> Dict:
    """
    Core DP analysis metrics calculation for filtered regions.

    Orchestrates DP analysis workflow: runs DP algorithm on each region,
    calculates uncertainties, propagates errors, and builds result metrics.

    Args:
        regions_dict: Dictionary of region_id -> variants list
        total_regions: Total number of regions in analysis
        uncertain_regions: Count of uncertain regions (N/A variants)
        variants_deterministic: Total count of variants for deterministic analysis
        unstable_threshold: P(≥1 MSI) threshold for calling region unstable
        msi_high_threshold: MSI score threshold for MSI-High classification
    """
    # Region counting
    regions_unstable_probabilistic = 0
    regions_with_variants = 0
    regions_unstable_expected = 0.0

    # Variants Counting
    variants_expected = 0.0

    # Uncertainty collection
    region_variances = []
    region_unstable_uncertainties = []

    # Fragile regions counting
    fragile_regions_count = 0

    for _, variants in regions_dict.items():
        if variants:
            regions_with_variants += 1

            # Core DP calculation
            distribution = run_msi_dp(variants)
            p_unstable = calculate_p_unstable(distribution)
            expected_value, variance = calculate_expected_and_variance(distribution)

            # Region classification
            if p_unstable >= unstable_threshold:
                regions_unstable_probabilistic += 1

            # Expected value accumulation
            regions_unstable_expected += p_unstable
            variants_expected += expected_value

            # Uncertainty collection
            region_variances.append(variance)
            bernoulli_variance = p_unstable * (1 - p_unstable)
            region_unstable_uncertainties.append(bernoulli_variance)

            # Fragile region counting
            if abs(p_unstable - unstable_threshold) <= 0.05:
                fragile_regions_count += 1


    # Core MSI score calculations
    msi_score_probabilistic = (regions_unstable_probabilistic / total_regions) * 100 if total_regions > 0 else 0.0
    msi_score_expected = (
        (regions_unstable_expected / total_regions) * 100 if total_regions > 0 else 0.0
    )
    msi_score_deterministic = (
        (regions_with_variants / total_regions) * 100 if total_regions > 0 else 0.0
    )

    # MSI status classifications
    msi_status_probabilistic = classify_msi_status(msi_score_probabilistic, msi_high_threshold)
    msi_status_expected = classify_msi_status(msi_score_expected, msi_high_threshold)
    msi_status_deterministic = classify_msi_status(
        msi_score_deterministic, msi_high_threshold
    )

    # Base uncertainty calculations
    uncertainty_unstable_expected = sum(region_unstable_uncertainties) ** 0.5
    uncertainty_variants_expected, _ = propagate_uncertainties(region_variances, total_regions)

    # Probabilistic method: Fragility analysis
    score_per_region_probabilistic = 100.0 / total_regions
    uncertainty_swing_probabilistic = fragile_regions_count * score_per_region_probabilistic
    uncertainty_range_probabilistic = [
        max(0.0, msi_score_probabilistic - uncertainty_swing_probabilistic),
        min(100.0, msi_score_probabilistic + uncertainty_swing_probabilistic)
    ]

    # MSI score uncertainties
    uncertainty_msi_score_expected = (
        (uncertainty_unstable_expected / total_regions) * 100 if total_regions > 0 else 0.0
    )

    # Impact and range calculations
    impact_msi_score_probabilistic_vs_deterministic = abs(msi_score_probabilistic - msi_score_deterministic)
    range_msi_score_probabilistic_vs_deterministic = [
        min(msi_score_probabilistic, msi_score_deterministic),
        max(msi_score_probabilistic, msi_score_deterministic),
    ]
    impact_unstable_expected_vs_deterministic = abs(regions_unstable_expected - regions_with_variants)
    range_unstable_expected_vs_deterministic = [
        min(regions_unstable_expected, regions_with_variants),
        max(regions_unstable_expected, regions_with_variants),
    ]
    range_msi_score_expected_uncertainty = [
        round(max(0.0, msi_score_expected - uncertainty_msi_score_expected), 2),
        round(min(100.0, msi_score_expected + uncertainty_msi_score_expected), 2)
    ]
    range_msi_score_overall_uncertainty = [
        round(min(uncertainty_range_probabilistic[0], range_msi_score_expected_uncertainty[0], msi_score_deterministic), 2),
        round(max(uncertainty_range_probabilistic[1], range_msi_score_expected_uncertainty[1], msi_score_deterministic), 2)
    ]
    range_unstable_expected_uncertainty = [
        round(max(0.0, regions_unstable_expected - uncertainty_unstable_expected), 1),
        round(regions_unstable_expected + uncertainty_unstable_expected, 1)
    ]
    range_unstable_probabilistic_fragility = [
        max(0, regions_unstable_probabilistic - fragile_regions_count),
        regions_unstable_probabilistic + fragile_regions_count
    ]
    range_unstable_overall_uncertainty = [
        min(range_unstable_expected_uncertainty[0], range_unstable_probabilistic_fragility[0], regions_with_variants),
        max(range_unstable_expected_uncertainty[1], range_unstable_probabilistic_fragility[1], regions_with_variants)
    ]
    range_variants_expected_uncertainty = [
        round(max(0.0, variants_expected - uncertainty_variants_expected), 2),
        round(variants_expected + uncertainty_variants_expected, 2)
    ]
    range_variants_overall_uncertainty = [
        round(min(range_variants_expected_uncertainty[0], variants_deterministic), 2),
        round(max(range_variants_expected_uncertainty[1], variants_deterministic), 2)
    ]

    # MSI variants calculations
    impact_variants_expected_vs_deterministic = abs(variants_expected - variants_deterministic)
    range_variants_expected_vs_deterministic = [
        min(variants_expected, variants_deterministic),
        max(variants_expected, variants_deterministic),
    ]

    # Region breakdown calculations
    regions_stable_deterministic = (
        total_regions - regions_with_variants - uncertain_regions
    )
    regions_stable_probabilistic = total_regions - regions_unstable_probabilistic - uncertain_regions

    return {
        #TODO: BETTER NESTED STRUCTURE
        # CORE MSI SCORE CALCULATIONS
        "msi_score_probabilistic": round(msi_score_probabilistic, 2),
        "msi_score_expected": round(msi_score_expected, 2),
        "msi_score_deterministic": round(msi_score_deterministic, 2),
        # MSI STATUS CLASSIFICATIONS
        "msi_status_probabilistic": msi_status_probabilistic,
        "msi_status_expected": msi_status_expected,
        "msi_status_deterministic": msi_status_deterministic,
        # BASE UNCERTAINTY CALCULATIONS
        "uncertainty_variants_expected": round(uncertainty_variants_expected, 2),
        "uncertainty_unstable_expected": round(uncertainty_unstable_expected, 2),
        # Probabilistic method uncertainty (fragility analysis)
        "fragile_regions_count": fragile_regions_count,
        "uncertainty_swing_msi_score_probabilistic": round(uncertainty_swing_probabilistic, 3),
        "uncertainty_range_msi_score_probabilistic": [round(x, 2) for x in uncertainty_range_probabilistic],
        "uncertainty_explanation_probabilistic": f"+/-{fragile_regions_count} regions within 0.05 of threshold",
        # MSI SCORE UNCERTAINTIES
        "uncertainty_msi_score_expected": round(uncertainty_msi_score_expected, 3),
        # IMPACT AND RANGE CALCULATIONS
        "impact_msi_score_probabilistic_vs_deterministic": round(impact_msi_score_probabilistic_vs_deterministic, 2),
        "range_msi_score_probabilistic_vs_deterministic": range_msi_score_probabilistic_vs_deterministic,
        "range_msi_score_expected_uncertainty": range_msi_score_expected_uncertainty,
        "range_msi_score_overall_uncertainty": range_msi_score_overall_uncertainty,
        "impact_unstable_expected_vs_deterministic": round(impact_unstable_expected_vs_deterministic, 1),
        "range_unstable_expected_vs_deterministic": range_unstable_expected_vs_deterministic,
        "range_unstable_expected_uncertainty": range_unstable_expected_uncertainty,
        "range_unstable_probabilistic_fragility": range_unstable_probabilistic_fragility,
        "range_unstable_overall_uncertainty": range_unstable_overall_uncertainty,
        # MSI VARIANTS CALCULATIONS
        "variants_expected": round(variants_expected, 2),
        "variants_deterministic": variants_deterministic,
        "impact_variants_expected_vs_deterministic": round(impact_variants_expected_vs_deterministic, 2),
        "range_variants_expected_vs_deterministic": range_variants_expected_vs_deterministic,
        "range_variants_expected_uncertainty": range_variants_expected_uncertainty,
        "range_variants_overall_uncertainty": range_variants_overall_uncertainty,
        # REGION BREAKDOWN CALCULATIONS
        # 1. Basic region counts
        "total_regions": total_regions,
        "regions_with_variants": regions_with_variants,
        "uncertain_regions": uncertain_regions,
        # 2. Method-specific region counts
        "regions_unstable_probabilistic": regions_unstable_probabilistic,
        "regions_unstable_expected": round(regions_unstable_expected, 2),
        "regions_unstable_deterministic": regions_with_variants,
        "regions_stable_probabilistic": regions_stable_probabilistic,
        "regions_stable_deterministic": regions_stable_deterministic,
        # Analysis Parameters
        "analysis_parameters": {
            "unstable_threshold": unstable_threshold,
            "msi_high_threshold": msi_high_threshold,
        }
    }


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

    for variant in variants:
        af_by_sample = variant["dp_data"]["af_by_sample"]
        sample_af = af_by_sample[sample_name]

        if af_threshold == -1:
            if sample_af == -1:
                filtered_variants.append(variant)
        elif af_threshold == -2:
            if sample_af == -2:
                filtered_variants.append(variant)
        else:
            if sample_af >= af_threshold:
                filtered_variants.append(variant)

    return filtered_variants


def run_af_evolution_analysis(
    dp_ready_data: Dict,
    total_ms_regions: int,
    unstable_threshold: float = 0.5,
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
    print("[AF-EVOLUTION] Starting AF evolution analysis")

    results = {}
    samples = dp_ready_data["samples"]
    af_thresholds = dp_ready_data["af_thresholds"]

    for sample_name in samples:
        print(f"[AF-EVOLUTION] Analyzing sample: {sample_name}")
        sample_results = {}

        for af_threshold in af_thresholds:
            af_data = dp_ready_data["af_data_by_sample"][sample_name][af_threshold]
            uncertain_regions = af_data["uncertain_regions"]
            total_variants = af_data["variant_count"]

            # Create filtered regions dictionary by sample and AF threshold
            # for calculating MSI metrics per AF threshold
            filtered_regions = {}
            for region_id, variants in dp_ready_data["regions"].items():
                filtered_variants = filter_variants_by_af_and_sample(
                    variants, af_threshold, sample_name
                )
                if filtered_variants:
                    filtered_regions[region_id] = filtered_variants

            metrics = calculate_msi_metrics_for_regions(
                filtered_regions,
                total_ms_regions,
                uncertain_regions,
                total_variants,
                unstable_threshold,
                msi_high_threshold,
            )

            sample_results[f"af_{af_threshold}"] = metrics

        results[sample_name] = sample_results

    print("[AF-EVOLUTION] AF evolution analysis complete")
    return results


def debug_af_filtering_discrepancy(
    dp_ready_data, total_ms_regions, unstable_threshold=0.5
):
    """Debug function to find the exact region causing the 1-region discrepancy"""

    print("\n" + "=" * 60)
    print("DEBUGGING AF FILTERING DISCREPANCY")
    print("=" * 60)

    # Regional analysis regions
    regional_unstable_regions = []
    for region_id, variants in dp_ready_data["regions"].items():
        if variants:
            distribution = run_msi_dp(variants)
            p_unstable = calculate_p_unstable(distribution)
            if p_unstable >= unstable_threshold:
                regional_unstable_regions.append(region_id)

    # AF >= 0.0 analysis regions
    af_0_unstable_regions = []
    sample_name = dp_ready_data["samples"][0]  # Use first sample

    for region_id, variants in dp_ready_data["regions"].items():
        filtered_variants = filter_variants_by_af_and_sample(variants, 0.0, sample_name)
        if filtered_variants:
            distribution = run_msi_dp(filtered_variants)
            p_unstable = calculate_p_unstable(distribution)
            if p_unstable >= unstable_threshold:
                af_0_unstable_regions.append(region_id)

    # Find the difference
    missing_regions = set(regional_unstable_regions) - set(af_0_unstable_regions)
    extra_regions = set(af_0_unstable_regions) - set(regional_unstable_regions)

    print(f"Regional unstable regions: {len(regional_unstable_regions)}")
    print(f"AF >= 0.0 unstable regions: {len(af_0_unstable_regions)}")
    print(f"Missing regions in AF >= 0.0: {len(missing_regions)}")
    print(f"Extra regions in AF >= 0.0: {len(extra_regions)}")

    # Analyze the missing region(s)
    for region_id in missing_regions:
        print(f"\nMISSING REGION: {region_id}")
        variants = dp_ready_data["regions"][region_id]
        print(f"  Total variants in region: {len(variants)}")

        for i, variant in enumerate(variants):
            af_data = variant["dp_data"]["af_by_sample"]
            print(f"  Variant {i}: AF = {af_data}")

        # Test filtering
        filtered = filter_variants_by_af_and_sample(variants, 0.0, sample_name)
        print(f"  Variants after AF >= 0.0 filtering: {len(filtered)}")

        if not filtered:
            print("  → This region was DROPPED because all variants have AF < 0.0!")
            # Check if all variants have AF = -1
            all_af_values = [
                variant["dp_data"]["af_by_sample"][sample_name] for variant in variants
            ]
            print(f"  → All AF values: {all_af_values}")

    print("=" * 60)
