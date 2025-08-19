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


def validate_probabilities(
    pp: Optional[float], pa: Optional[float], art: Optional[float]
) -> Tuple[bool, List[str]]:
    """
    Validate probability values for present, absent, and artifact states.

    Args:
        pp: Probability present (0-1 or None)
        pa: Probability absent (0-1 or None)
        art: Probability artifact (0-1 or None)

    Returns:
        Tuple of (is_valid, error_list)

    Note:
        Uses 0.005 tolerance for floating-point precision in PHRED-to-probability conversion.
    """
    errors = []

    for name, value in [("present", pp), ("absent", pa), ("artifact", art)]:
        if value is not None:
            if not isinstance(value, (int, float)):
                errors.append(f"invalid_type_{name}")
                continue
            if value < 0 or value > 1:
                errors.append(f"invalid_range_{name}_{value}")

    known_values = [x for x in [pp, pa, art] if x is not None]
    if known_values:
        known_sum = sum(known_values)
        if known_sum > 1.005:
            errors.append(f"sum_exceeds_1.0_{known_sum:.6f}")

    return len(errors) == 0, errors


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
def apply_4_step_imputation_to_variant(variant: Dict) -> None:
    """
    Apply 4-step probabilistic imputation to variant with complete audit trail.

    Implements formula for MSI analysis: p_present + p_absent = 1.0
    where p_absent = prob_absent + prob_artifact (consolidation step).

    IMPUTATION STRATEGY:
    1. All present (0 missing): Direct consolidation
    2. One missing: Derive from constraint (sum = 1.0)
    3. Two missing: Proportional split using 2:1:1 ratio (pp:pa:art)
    4. All missing: Uniform distribution (present=0.5, absent=0.5)

    Args:
        variant (Dict): Variant dictionary containing prob_present, prob_absent,
                        prob_artifact, and sample_afs fields

    Returns:
        None: Modifies variant in-place by adding dp_data field with:
            - p_present: Final present probability (0-1)
            - p_absent: Final absent probability (0-1)
            - af_by_sample: Validated AF values per sample (-2=invalid, -1=missing, 0-1=valid)
            - audit_trail: Complete imputation methodology tracking

    Side Effects:
        - Adds variant["dp_data"] field with imputation results
        - Validates and normalizes all AF values with error tracking
        - Creates detailed audit trail for reproducibility
    """

    pp = variant.get("prob_present")
    pa = variant.get("prob_absent")
    art = variant.get("prob_artifact")

    audit_trail = {
        "missing_fields": [],
        "imputation_method": None,
        "calculation_steps": None,
        "validation_issues": [],
        "fallback_applied": False,
        "af_processing": {
            "missing_af_samples": [],
            "af_conversion": {},
            "af_errors": [],
        },
    }

    is_valid, prob_errors = validate_probabilities(pp, pa, art)
    if not is_valid:
        p_present = 0.5
        p_absent = 0.5
        audit_trail["imputation_method"] = "validation_fallback_uniform"
        audit_trail["calculation_steps"] = (
            f"validation failed: {prob_errors}, "
            f"applied uniform p_present=0.5, p_absent=0.5"
        )
        audit_trail["validation_issues"] = prob_errors
        audit_trail["fallback_applied"] = True
    else:
        missing_fields = []
        if pp is None:
            missing_fields.append("prob_present")
        if pa is None:
            missing_fields.append("prob_absent")
        if art is None:
            missing_fields.append("prob_artifact")

        audit_trail["missing_fields"] = missing_fields
        missing_count = len(missing_fields)

        # SECTION: All probabilities available - direct consolidation
        if missing_count == 0:
            p_present = pp
            p_absent = pa + art
            audit_trail["imputation_method"] = "consolidation_direct"
            audit_trail["calculation_steps"] = (
                f"p_absent = prob_absent + prob_artifact = {pa} + {art} = {p_absent}, "
                f"p_present = {p_present}"
            )

        # SECTION: One missing - derive from constraint (sum = 1.0)
        elif missing_count == 1:
            if pp is None:
                derived = 1.0 - pa - art
                p_present = derived
                p_absent = pa + art
                audit_trail["imputation_method"] = "constraint_derived_present"
                audit_trail["calculation_steps"] = (
                    f"derived prob_present = 1.0 - {pa} - {art} = {derived}, "
                    f"p_absent = {pa} + {art} = {p_absent}"
                )

            elif pa is None:
                derived = 1.0 - pp - art
                p_present = pp
                p_absent = derived + art
                audit_trail["imputation_method"] = "constraint_derived_absent"
                audit_trail["calculation_steps"] = (
                    f"derived prob_absent = 1.0 - {pp} - {art} = {derived}, "
                    f"p_absent = {derived} + {art} = {p_absent}"
                )

            elif art is None:
                derived = 1.0 - pp - pa
                p_present = pp
                p_absent = pa + derived
                audit_trail["imputation_method"] = "constraint_derived_artifact"
                audit_trail["calculation_steps"] = (
                    f"derived prob_artifact = 1.0 - {pp} - {pa} = {derived}, "
                    f"p_absent = {pa} + {derived} = {p_absent}"
                )

        # SECTION: All missing - uniform distribution
        elif missing_count == 3:
            p_present = 0.5
            p_absent = 0.5
            audit_trail["imputation_method"] = "uniform_all_missing"
            audit_trail["calculation_steps"] = (
                "all probabilities missing, applied uniform: present=0.5, absent=0.25, "
                "artifact=0.25, p_absent=0.5"
            )

        # SECTION: Two missing - proportional split
        # f1: (pa + art + pp = 1.0)
        # f2: (pa + art = p_absent)
        # Assumed amputations: pp=0.5, pa=0.25, art=0.25
        # Ratio is: pp:pa:art = 2:1:1
        else:
            known_value = pp if pp is not None else (pa if pa is not None else art)
            remaining = 1.0 - known_value

            if pp is None and pa is None:
                imputed_present = remaining * (2 / 3)
                imputed_absent = remaining * (1 / 3)
                p_present = imputed_present
                p_absent = imputed_absent + art
                audit_trail["imputation_method"] = "proportional_split_present_absent"
                audit_trail["calculation_steps"] = (
                    f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, "
                    f"absent = {remaining} * (1/3) = {imputed_absent}, "
                    f"p_absent = {imputed_absent} + {art} = {p_absent}"
                )

            elif pp is None and art is None:
                imputed_present = remaining * (2 / 3)
                imputed_artifact = remaining * (1 / 3)
                p_present = imputed_present
                p_absent = pa + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_present_artifact"
                audit_trail["calculation_steps"] = (
                    f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, "
                    f"artifact = {remaining} * (1/3) = {imputed_artifact}, "
                    f"p_absent = {pa} + {imputed_artifact} = {p_absent}"
                )

            elif pa is None and art is None:
                imputed_absent = remaining * 0.5
                imputed_artifact = remaining * 0.5
                p_present = pp
                p_absent = imputed_absent + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_absent_artifact"
                audit_trail["calculation_steps"] = (
                    f"remaining = {remaining}, absent = {remaining} * 0.5 = {imputed_absent}, "
                    f"artifact = {remaining} * 0.5 = {imputed_artifact}, "
                    f"p_absent = {imputed_absent} + {imputed_artifact} = {p_absent}"
                )

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

    #TODO: PERFORMANCE OPTIMIZATION CONSIDERATION
    # Current: Separate loops of same kind (data prep + imputation)
    # Future: Consider merging into single func with AF logic
    # Trade-off: Performance gain vs code clarity/maintainability
    # State: Profile first, optimize only if bottleneck identified
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
                # Apply 4-step imputation
                apply_4_step_imputation_to_variant(variant)

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


def calculate_variance(distribution: List[float]) -> float:
    """
    Calculate variance from DP probability distribution.

    Mathematical foundation: Var(X) = Σ((i - E[X])² × P(X = i)) for i = 0 to n
    This represents the variance in the number of MSI variants.

    Args:
        distribution (List[float]): Probability distribution from run_msi_dp()
            where distribution[i] = P(exactly i MSI variants)

    Returns:
        float: Variance in number of MSI variants

    Note:
        Used for DP uncertainty propagation in expected unstable regions and
        expected MSI variants. Distribution validation handled by calculate_expected_value.
    """
    expected = calculate_expected_value(distribution)
    variance = sum((i - expected) ** 2 * prob for i, prob in enumerate(distribution))
    return variance


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
    total_variants: int,
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
        total_variants: Total count of variants for deterministic_msi_variants
        unstable_threshold: P(≥1 MSI) threshold for calling region unstable
        msi_high_threshold: MSI score threshold for MSI-High classification
    """
    # Region counting
    unstable_count = 0
    regions_with_variants = 0
    expected_unstable_regions = 0.0

    # Variants Counting
    expected_msi_variants = 0.0

    # Uncertainty collection
    region_variances = []
    region_unstable_uncertainties = []

    # uncertain_regions_count = dp_ready_data.get("total_uncertain_regions", 0)

    for region_id, variants in regions_dict.items():
        if variants:
            regions_with_variants += 1

            # Core DP calculation
            distribution = run_msi_dp(variants)
            p_unstable = calculate_p_unstable(distribution)
            variance = calculate_variance(distribution)

            # Region classification
            if p_unstable > unstable_threshold:
                unstable_count += 1

            # Expected value accumulation
            expected_unstable_regions += p_unstable
            expected_msi_variants += calculate_expected_value(distribution)

            # Uncertainty collection
            region_variances.append(variance)
            bernoulli_variance = p_unstable * (1 - p_unstable)
            region_unstable_uncertainties.append(bernoulli_variance)

    # Core MSI score calculations
    msi_score = (unstable_count / total_regions) * 100 if total_regions > 0 else 0.0
    expected_msi_score = (
        (expected_unstable_regions / total_regions) * 100 if total_regions > 0 else 0.0
    )
    deterministic_msi_score = (
        (regions_with_variants / total_regions) * 100 if total_regions > 0 else 0.0
    )

    # MSI status classifications
    msi_status = classify_msi_status(msi_score, msi_high_threshold)
    expected_msi_status = classify_msi_status(expected_msi_score, msi_high_threshold)
    deterministic_msi_status = classify_msi_status(
        deterministic_msi_score, msi_high_threshold
    )

    # Base uncertainty calculations
    unstable_regions_uncertainty = sum(region_unstable_uncertainties) ** 0.5
    absolute_uncertainty, overall_uncertainty = propagate_uncertainties(
        region_variances, total_regions
    )
    expected_variants_uncertainty = absolute_uncertainty

    # MSI score uncertainties
    msi_score_statistical_uncertainty = (
        (unstable_regions_uncertainty / total_regions) * 100
        if total_regions > 0
        else 0.0
    )
    expected_msi_score_uncertainty = (
        (absolute_uncertainty / total_regions) * 100 if total_regions > 0 else 0.0
    )

    # Impact and range calculations
    analysis_impact = abs(msi_score - deterministic_msi_score)
    msi_score_range = [
        min(msi_score, deterministic_msi_score),
        max(msi_score, deterministic_msi_score),
    ]
    expected_unstable_impact = abs(expected_unstable_regions - regions_with_variants)
    expected_unstable_range = [
        min(expected_unstable_regions, regions_with_variants),
        max(expected_unstable_regions, regions_with_variants),
    ]

    # MSI variants calculations
    deterministic_msi_variants = total_variants
    msi_variants_impact = abs(expected_msi_variants - deterministic_msi_variants)
    msi_variants_range = [
        min(expected_msi_variants, deterministic_msi_variants),
        max(expected_msi_variants, deterministic_msi_variants),
    ]

    # Region breakdown calculations
    deterministic_stable_regions = (
        total_regions - regions_with_variants - uncertain_regions
    )
    probabilistic_stable_regions = total_regions - unstable_count - uncertain_regions
    deterministic_unstable_regions = regions_with_variants

    return {
        # CORE MSI SCORE CALCULATIONS
        "msi_score_probabilistic": round(msi_score, 2),
        "msi_score_expected": round(expected_msi_score, 2),
        "msi_score_deterministic": round(deterministic_msi_score, 2),
        # MSI STATUS CLASSIFICATIONS
        "msi_status_probabilistic": msi_status,
        "msi_status_expected": expected_msi_status,
        "msi_status_deterministic": deterministic_msi_status,
        # BASE UNCERTAINTY CALCULATIONS
        "uncertainty_variants_expected": round(expected_variants_uncertainty, 2),
        "uncertainty_unstable_expected": round(absolute_uncertainty, 2),
        # MSI SCORE UNCERTAINTIES
        "uncertainty_msi_score_probabilistic": round(
            msi_score_statistical_uncertainty, 3
        ),
        "uncertainty_msi_score_expected": round(expected_msi_score_uncertainty, 3),
        # IMPACT AND RANGE CALCULATIONS
        "impact_msi_score_probabilistic_vs_deterministic": round(analysis_impact, 2),
        "range_msi_score_probabilistic_vs_deterministic": msi_score_range,
        "impact_unstable_expected_vs_deterministic": round(expected_unstable_impact, 1),
        "range_unstable_expected_vs_deterministic": expected_unstable_range,
        # MSI VARIANTS CALCULATIONS
        "variants_expected": round(expected_msi_variants, 2),
        "variants_deterministic": deterministic_msi_variants,
        "impact_variants_expected_vs_deterministic": round(msi_variants_impact, 2),
        "range_variants_expected_vs_deterministic": msi_variants_range,
        # REGION BREAKDOWN CALCULATIONS
        # 1. Basic region counts
        "total_regions": total_regions,
        "regions_with_variants": regions_with_variants,
        "regions_without_variants": total_regions - regions_with_variants,
        "uncertain_regions": uncertain_regions,
        # 2. Method-specific region counts
        "regions_unstable_probabilistic": unstable_count,
        "regions_unstable_expected": round(expected_unstable_regions, 2),
        "regions_unstable_deterministic": deterministic_unstable_regions,
        "regions_stable_probabilistic": probabilistic_stable_regions,
        "regions_stable_deterministic": deterministic_stable_regions,
    }


def run_regional_msi_analysis(
    dp_ready_data: Dict,
    total_ms_regions_from_bed: int,
    unstable_threshold: float = 0.5,
    msi_high_threshold: float = 3.5,
) -> Dict:
    """
    Run genome-wide regional MSI analysis using dynamic programming approach.

    Wrapper function that calls calculate_msi_metrics_for_regions() with
    genome-wide data from prepare_variants_for_dp().

    Args:
        dp_ready_data: Prepared data from prepare_variants_for_dp()
        total_ms_regions_from_bed: Total microsatellite regions from BED file
        unstable_threshold: P(≥1 MSI) threshold for region classification (default: 0.5)
        msi_high_threshold: MSI score threshold for MSI-High classification (default: 3.5)

    Returns:
        Dict: Complete MSI analysis results from calculate_msi_metrics_for_regions()
    """
    print("[REGIONAL-DP] Starting regional MSI analysis")

    return calculate_msi_metrics_for_regions(
        dp_ready_data["regions"],
        total_ms_regions_from_bed,
        dp_ready_data.get("total_uncertain_regions", 0),
        dp_ready_data["total_variants"],
        unstable_threshold,
        msi_high_threshold,
    )


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

        af_scores = [sample_results[f"af_{af}"]["msi_score"] for af in af_thresholds]
        print(
            f"[AF-EVOLUTION] {sample_name} evolution: {af_thresholds[0]}→{af_thresholds[-1]}: {af_scores[0]:.1f}%→{af_scores[-1]:.1f}%"
        )

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
            if p_unstable > unstable_threshold:
                regional_unstable_regions.append(region_id)

    # AF >= 0.0 analysis regions
    af_0_unstable_regions = []
    sample_name = dp_ready_data["samples"][0]  # Use first sample

    for region_id, variants in dp_ready_data["regions"].items():
        filtered_variants = filter_variants_by_af_and_sample(variants, 0.0, sample_name)
        if filtered_variants:
            distribution = run_msi_dp(filtered_variants)
            p_unstable = calculate_p_unstable(distribution)
            if p_unstable > unstable_threshold:
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
