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


def apply_4_step_imputation_to_variant(variant: Dict) -> None:
    """Apply complete 4-step imputation with Johannes formula and full audit trail."""

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
            f"validation failed: {prob_errors}, applied uniform p_present=0.5, p_absent=0.5"
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
                f"p_absent = prob_absent + prob_artifact = {pa} + {art} = {p_absent}, p_present = {p_present}"
            )

        # SECTION: One missing - derive from constraint (sum = 1.0)  
        elif missing_count == 1:
            if pp is None:
                derived = 1.0 - pa - art
                p_present = derived
                p_absent = pa + art
                audit_trail["imputation_method"] = "constraint_derived_present"
                audit_trail["calculation_steps"] = (
                    f"derived prob_present = 1.0 - {pa} - {art} = {derived}, p_absent = {pa} + {art} = {p_absent}"
                )

            elif pa is None:
                derived = 1.0 - pp - art
                p_present = pp
                p_absent = derived + art
                audit_trail["imputation_method"] = "constraint_derived_absent"
                audit_trail["calculation_steps"] = (
                    f"derived prob_absent = 1.0 - {pp} - {art} = {derived}, p_absent = {derived} + {art} = {p_absent}"
                )

            elif art is None:
                derived = 1.0 - pp - pa
                p_present = pp
                p_absent = pa + derived
                audit_trail["imputation_method"] = "constraint_derived_artifact"
                audit_trail["calculation_steps"] = (
                    f"derived prob_artifact = 1.0 - {pp} - {pa} = {derived}, p_absent = {pa} + {derived} = {p_absent}"
                )

        # SECTION: All missing - uniform distribution
        elif missing_count == 3:
            p_present = 0.5
            p_absent = 0.5
            audit_trail["imputation_method"] = "uniform_all_missing"
            audit_trail["calculation_steps"] = (
                "all probabilities missing, applied uniform: present=0.5, absent=0.25, artifact=0.25, p_absent=0.5"
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
                    f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, absent = {remaining} * (1/3) = {imputed_absent}, p_absent = {imputed_absent} + {art} = {p_absent}"
                )

            elif pp is None and art is None:
                imputed_present = remaining * (2 / 3)
                imputed_artifact = remaining * (1 / 3)
                p_present = imputed_present
                p_absent = pa + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_present_artifact"
                audit_trail["calculation_steps"] = (
                    f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, artifact = {remaining} * (1/3) = {imputed_artifact}, p_absent = {pa} + {imputed_artifact} = {p_absent}"
                )

            elif pa is None and art is None:
                imputed_absent = remaining * 0.5
                imputed_artifact = remaining * 0.5
                p_present = pp
                p_absent = imputed_absent + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_absent_artifact"
                audit_trail["calculation_steps"] = (
                    f"remaining = {remaining}, absent = {remaining} * 0.5 = {imputed_absent}, artifact = {remaining} * 0.5 = {imputed_artifact}, p_absent = {imputed_absent} + {imputed_artifact} = {p_absent}"
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


def extract_af_data_for_variant(variant, sample_list):
    """
    Extract and normalize AF data without full DP processing
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


def prepare_variants_for_dp(results, vcf_file_path, imputation_method="uniform"):
    """Data Preparation Function for Dynamic Programming Analysis."""

    print("[DP-ENTRY] Starting data preparation for dynamic programming analysis")
    print(f"[DP-ENTRY] Processing {len(results)} regions")

    vcf = pysam.VariantFile(vcf_file_path)
    sample_list = list(vcf.header.samples)
    vcf.close()
    print(f"[DP-PREP] Found samples in VCF: {sample_list}")

    af_thresholds = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1, -2]

    total_uncertain_regions = 0

    uncertain_by_af_by_sample = {}
    for sample_name in sample_list:
        uncertain_by_af_by_sample[sample_name] = {
            af_threshold: 0 for af_threshold in af_thresholds
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
                        sample_af = af_data[sample_name]
                        if (
                            (af_threshold == -1 and sample_af == -1)
                            or (af_threshold == -2 and sample_af == -2)
                            or (af_threshold >= 0 and sample_af >= af_threshold)
                        ):
                            region_has_variants_at_af_for_sample = True
                            break

                    if region_has_variants_at_af_for_sample:
                        uncertain_by_af_by_sample[sample_name][af_threshold] += 1

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
        "uncertain_by_af_by_sample": uncertain_by_af_by_sample,
        "total_uncertain_regions": total_uncertain_regions,
        # "entry_summary": {
        #     "input_regions": len(results),
        #     "regions_with_variants": len(region_variants),
        #     "total_variants": total_variants,
        #     "step_2_complete": True
        # }
    }

    print(
        f"[DP-PREP] Processed {total_variants} variants in {len(region_variants)} regions"
    )
    print(f"[DP-PREP] Found {total_uncertain_regions} uncertain regions")
    print(f"[DP-PREP] AF thresholds: {af_thresholds}")
    print("[DP-PREP] Data Preparation for DP complete")

    return dp_ready


def run_msi_dp(variants_with_probabilities):
    """
    Matrix m: n variants (columns) × n+1 possible MSI counts (rows)
    m[i,k] = P(exactly i MSI variants using first k variants)

    Args:
        variants_with_probabilities: List of variants with dp_data containing
            p_present and p_absent probabilities

    Returns:
        List: [P(0 MSI), P(1 MSI), P(2 MSI), ..., P(n MSI)]
    """
    n = len(variants_with_probabilities)

    if n == 0:
        return [1.0]

    prev_col = [0.0] * (n + 1)

    p_absent_0 = variants_with_probabilities[0]["dp_data"]["p_absent"]
    p_present_0 = variants_with_probabilities[0]["dp_data"]["p_present"]

    prev_col[0] = p_absent_0
    prev_col[1] = p_present_0

    for k in range(1, n):
        p_absent_k = variants_with_probabilities[k]["dp_data"]["p_absent"]
        p_present_k = variants_with_probabilities[k]["dp_data"]["p_present"]

        curr_col = [0.0] * (n + 1)

        # Base case: P(0 MSI) = previous P(0 MSI) × variant absent
        curr_col[0] = prev_col[0] * p_absent_k

        # Recurrence: Formula for rows 1 to k+1
        # m[i,k] = m[i,k-1] × p_absent + m[i-1,k-1] × p_present
        for i in range(1, k + 2):
            curr_col[i] = prev_col[i] * p_absent_k + prev_col[i - 1] * p_present_k

        prev_col = curr_col

    return prev_col


def calculate_expected_value(distribution):
    """
    Calculate expected number of MSI variants from probability distribution.
    """
    return sum(i * prob for i, prob in enumerate(distribution))


def calculate_std_dev(distribution):
    """
    Calculate standard deviation (uncertainty) from probability distribution.
    """
    expected = calculate_expected_value(distribution)
    variance = sum((i - expected) ** 2 * prob for i, prob in enumerate(distribution))
    return variance**0.5


def calculate_p_unstable(distribution):
    """Calculate probability of instability: P(≥1 MSI variant)."""
    return 1.0 - distribution[0]


def propagate_uncertainties(uncertainties, total_regions):
    """
    Standard statistical uncertainty propagation for MSI score.

    MATHEMATICAL FOUNDATION:
    For independent random variables X₁, X₂, ..., Xₙ:
    - Var(X₁ + X₂ + ... + Xₙ) = Var(X₁) + Var(X₂) + ... + Var(Xₙ)
    - σ(sum) = √(σ₁² + σ₂² + ... + σₙ²)

    APPLICATION TO MSI:
    - Each region contributes: P(unstable) ± uncertainty to total unstable count
    - Total unstable count = sum of contributions ± √(sum of uncertainty²)
    - MSI score = (total unstable / total regions) × 100
    - MSI uncertainty = (uncertainty in count / total regions) × 100

    ERROR PROPAGATION FOR DIVISION:
    If Y = X/N (constant), then σ_Y = σ_X/N

    Args:
        uncertainties: List of standard deviations from individual regions
        total_regions: Total MS regions in genome (denominator for MSI score)

    Returns:
        float: Uncertainty in MSI score (percentage points)
    """
    if not uncertainties:
        return 0.0, 0.0

    # Sum of variances
    total_variance = sum(u**2 for u in uncertainties)

    # Standard deviation of total unstable count
    uncertainty_in_count = total_variance**0.5

    # Convert to MSI score uncertainty (error propagation for division)
    # MSI uncertainty = σ_X / N
    msi_uncertainty = (uncertainty_in_count / total_regions) * 100

    return uncertainty_in_count, msi_uncertainty


def classify_msi_status(msi_score, msi_high_threshold=3.5):
    """Classify MSI status using MSIsensor standard."""
    return "MSI-High" if msi_score >= msi_high_threshold else "MSS"


def run_regional_msi_analysis(
    dp_ready_data,
    total_ms_regions_from_bed,
    unstable_threshold=0.5,
    msi_high_threshold=3.5,
):
    """
    Regional MSI analysis using DP method.

    Args:
        dp_ready_data: Output from prepare_variants_for_dp()
        unstable_threshold: P(≥1 MSI) threshold for unstable classification
        msi_high_threshold: MSI score threshold for MSI-High classification
        total_ms_regions_from_bed: Total MS regions from BED file

    Returns:
        dict: Single MSI analysis result
    """

    print("[REGIONAL-DP] Starting regional MSI analysis")

    unstable_count = 0
    region_uncertainties = []
    regions_with_variants = 0
    uncertain_regions_count = dp_ready_data.get("total_uncertain_regions", 0)

    expected_unstable_regions = 0.0
    expected_msi_variants = 0.0
    region_unstable_uncertainties = []

    for region_id, variants in dp_ready_data["regions"].items():
        if variants:
            regions_with_variants += 1

            distribution = run_msi_dp(variants)
            p_unstable = calculate_p_unstable(distribution)
            uncertainty = calculate_std_dev(distribution)

            region_uncertainties.append(uncertainty)

            if p_unstable > unstable_threshold:
                unstable_count += 1

            expected_unstable_regions += p_unstable
            expected_msi_variants += calculate_expected_value(distribution)

            bernoulli_variance = p_unstable * (1 - p_unstable)
            region_unstable_uncertainties.append(bernoulli_variance)

    # Unstable regions uncertainty
    unstable_regions_uncertainty = sum(region_unstable_uncertainties) ** 0.5

    msi_score = (
        (unstable_count / total_ms_regions_from_bed) * 100
        if total_ms_regions_from_bed > 0
        else 0.0
    )
    msi_score_statistical_uncertainty = (
        (unstable_regions_uncertainty / total_ms_regions_from_bed) * 100
        if total_ms_regions_from_bed > 0
        else 0.0
    )
    absolute_uncertainty, overall_uncertainty = propagate_uncertainties(
        region_uncertainties, total_ms_regions_from_bed
    )
    msi_status = classify_msi_status(msi_score, msi_high_threshold)
    expected_variants_uncertainty = absolute_uncertainty

    # Calculate Deterministic MSI score and impact
    deterministic_msi_score = (
        (regions_with_variants / total_ms_regions_from_bed) * 100
        if total_ms_regions_from_bed > 0
        else 0.0
    )
    analysis_impact = abs(msi_score - deterministic_msi_score)
    msi_score_range = [
        min(msi_score, deterministic_msi_score),
        max(msi_score, deterministic_msi_score),
    ]
    deterministic_msi_status = classify_msi_status(
        deterministic_msi_score, msi_high_threshold
    )

    expected_unstable_impact = abs(expected_unstable_regions - regions_with_variants)
    expected_unstable_range = [
        min(expected_unstable_regions, regions_with_variants),
        max(expected_unstable_regions, regions_with_variants),
    ]

    # Deterministic MSI Variants
    deterministic_msi_variants = dp_ready_data["total_variants"]
    msi_variants_impact = abs(expected_msi_variants - deterministic_msi_variants)
    msi_variants_range = [
        min(expected_msi_variants, deterministic_msi_variants),
        max(expected_msi_variants, deterministic_msi_variants),
    ]

    deterministic_stable_regions = (
        total_ms_regions_from_bed - regions_with_variants - uncertain_regions_count
    )
    probabilistic_stable_regions = (
        total_ms_regions_from_bed - unstable_count - uncertain_regions_count
    )
    probabilistic_unstable_regions = unstable_count
    deterministic_unstable_regions = regions_with_variants

    print(
        f"[REGIONAL-DP] Results: {unstable_count}/{total_ms_regions_from_bed} unstable regions"
    )
    print(
        f"[REGIONAL-DP] MSI Score Range: {msi_score_range[0]:.2f}% - {msi_score_range[1]:.2f}% (methodological uncertainty)"
    )
    print(
        f"[REGIONAL-DP] Probabilistic MSI Score: {msi_score:.2f}% (P(≥1 MSI) > {unstable_threshold})"
    )
    print(
        f"[REGIONAL-DP] MSI Score Statistical Uncertainty: ± {msi_score_statistical_uncertainty:.3f}%"
    )
    print(
        f"[REGIONAL-DP] Deterministic MSI Score: {deterministic_msi_score:.2f}% (regions with ≥1 variant)"
    )
    print(f"[REGIONAL-DP] Analysis Impact: {analysis_impact:.2f} percentage points")
    print(f"[REGIONAL-DP] MSI Status (Probabilistic): {msi_status}")
    print(f"[REGIONAL-DP] MSI Status (Deterministic): {deterministic_msi_status}")
    print(f"[REGIONAL-DP] Total MS regions in BED: {total_ms_regions_from_bed}")
    print(f"[REGIONAL-DP] Regions with variants: {regions_with_variants}")
    print(f"[REGIONAL-DP] Uncertain regions: {uncertain_regions_count}")
    print(f"[REGIONAL-DP] Deterministic stable regions: {deterministic_stable_regions}")
    print(
        f"[REGIONAL-DP] Probabilistic unstable regions: {probabilistic_unstable_regions}"
    )
    print(f"[REGIONAL-DP] Probabilistic stable regions: {probabilistic_stable_regions}")
    print(
        f"[REGIONAL-DP] Deterministic unstable regions: {deterministic_unstable_regions}"
    )
    print(
        f"[REGIONAL-DP] Expected Unstable Regions Range: {expected_unstable_range[0]:.2f} - {expected_unstable_range[1]:.2f} (methodological uncertainty)"
    )
    print(
        f"[REGIONAL-DP] Probabilistic: {expected_unstable_regions:.2f}, Deterministic: {regions_with_variants}, Impact: {expected_unstable_impact:.2f}"
    )
    print(
        f"[REGIONAL-DP] Expected MSI Variants Range: {msi_variants_range[0]:.1f} - {msi_variants_range[1]:.1f} (methodological uncertainty)"
    )
    print(
        f"[REGIONAL-DP] Probabilistic: {expected_msi_variants:.2f} ± {expected_variants_uncertainty:.2f}, Deterministic: {deterministic_msi_variants}, Impact: {msi_variants_impact:.1f}"
    )
    print(
        f"[REGIONAL-DP] Expected Unstable Regions: {expected_unstable_regions:.2f} ± {unstable_regions_uncertainty:.2f}, Deterministic: {regions_with_variants}, Impact: {expected_unstable_impact:.2f}"
    )

    return {
        "msi_score": round(msi_score, 2),
        "msi_score_statistical_uncertainty": round(
            msi_score_statistical_uncertainty, 3
        ),
        # "msi_uncertainty": round(overall_uncertainty, 2),
        "msi_status": msi_status,
        "msi_score_range": msi_score_range,
        "msi_score_deterministic": round(deterministic_msi_score, 2),
        "analysis_impact": round(analysis_impact, 2),
        "msi_status_deterministic": deterministic_msi_status,
        "unstable_regions": unstable_count,
        "total_regions": total_ms_regions_from_bed,
        "regions_with_variants": regions_with_variants,
        "regions_without_variants": total_ms_regions_from_bed - regions_with_variants,
        "expected_unstable_regions": round(expected_unstable_regions, 2),
        "expected_unstable_range": expected_unstable_range,
        "expected_unstable_impact": round(expected_unstable_impact, 1),
        "expected_variants_uncertainty": round(expected_variants_uncertainty, 2),
        "expected_unstable_uncertainty": round(unstable_regions_uncertainty, 2),
        "expected_msi_variants": round(expected_msi_variants, 2),
        "deterministic_msi_variants": deterministic_msi_variants,
        "msi_variants_range": msi_variants_range,
        "msi_variants_impact": round(msi_variants_impact, 2),
        "uncertain_regions": uncertain_regions_count,
        "deterministic_stable_regions": deterministic_stable_regions,
        "probabilistic_unstable_regions": probabilistic_unstable_regions,
        "probabilistic_stable_regions": probabilistic_stable_regions,
        "deterministic_unstable_regions": deterministic_unstable_regions,
    }


####


def filter_variants_by_af_and_sample(variants, af_threshold, sample_name):
    """
    Filter variants by AF threshold for a specific sample.

    Args:
        variants: List of variants in a region
        af_threshold: Minimum AF threshold
        sample_name: Sample to filter for

    Returns:
        List of variants that meet AF threshold for this sample
    """
    filtered = []

    for variant in variants:
        af_by_sample = variant["dp_data"]["af_by_sample"]
        sample_af = af_by_sample[sample_name]

        if af_threshold == -1:
            if sample_af == -1:
                filtered.append(variant)
        elif af_threshold == -2:
            if sample_af == -2:
                filtered.append(variant)
        else:
            if sample_af >= af_threshold:
                filtered.append(variant)

    return filtered


def run_af_evolution_analysis(
    dp_ready_data, total_ms_regions, unstable_threshold=0.5, msi_high_threshold=3.5
):
    """
    AF-based MSI evolution analysis showing subclonal timeline.

    Shows how MSI score changes across AF thresholds, representing
    tumor evolution from early (high AF) to late (low AF) mutations.

    Args:
        dp_ready_data: Output from prepare_variants_for_dp()
        total_ms_regions: Total MS regions from BED file
        unstable_threshold: P(≥1 MSI) threshold for calling region unstable

    Returns:
        dict: AF evolution results grouped by sample
        {
            "sample1": {
                "af_1.0": {"msi_score": 2.1, "msi_uncertainty": 0.3},
                "af_0.8": {"msi_score": 4.7, "msi_uncertainty": 0.5},
                ...
            }
        }
    """

    print("[AF-EVOLUTION] Starting AF evolution analysis")

    results = {}
    samples = dp_ready_data["samples"]
    af_thresholds = dp_ready_data["af_thresholds"]

    for sample_name in samples:
        print(f"[AF-EVOLUTION] Analyzing sample: {sample_name}")
        sample_results = {}

        sample_uncertain_by_af = dp_ready_data.get("uncertain_by_af_by_sample", {}).get(
            sample_name, {}
        )

        for af_threshold in af_thresholds:
            unstable_count = 0
            region_uncertainties = []
            regions_analyzed = 0
            expected_unstable_regions = 0.0
            expected_msi_variants = 0.0
            af_unstable_uncertainties = []

            for region_id, variants in dp_ready_data["regions"].items():
                filtered_variants = filter_variants_by_af_and_sample(
                    variants, af_threshold, sample_name
                )

                if filtered_variants:
                    regions_analyzed += 1

                    distribution = run_msi_dp(filtered_variants)
                    p_unstable = calculate_p_unstable(distribution)
                    uncertainty = calculate_std_dev(distribution)

                    region_uncertainties.append(uncertainty)

                    if p_unstable > unstable_threshold:
                        unstable_count += 1

                    expected_unstable_regions += p_unstable
                    expected_msi_variants += calculate_expected_value(distribution)

                    bernoulli_variance = p_unstable * (1 - p_unstable)
                    af_unstable_uncertainties.append(bernoulli_variance)

            msi_score = (
                (unstable_count / total_ms_regions) * 100
                if total_ms_regions > 0
                else 0.0
            )
            absolute_uncertainty, msi_uncertainty = propagate_uncertainties(
                region_uncertainties, total_ms_regions
            )
            expected_variants_uncertainty = absolute_uncertainty

            af_unstable_regions_uncertainty = sum(af_unstable_uncertainties) ** 0.5
            af_msi_score_statistical_uncertainty = (
                (af_unstable_regions_uncertainty / total_ms_regions) * 100
                if total_ms_regions > 0
                else 0.0
            )

            # Deterministic MSI score
            deterministic_msi_score = (regions_analyzed / total_ms_regions) * 100

            analysis_impact = abs(msi_score - deterministic_msi_score)
            msi_score_range = [
                min(msi_score, deterministic_msi_score),
                max(msi_score, deterministic_msi_score),
            ]

            expected_unstable_impact = abs(expected_unstable_regions - regions_analyzed)
            expected_unstable_range = [
                min(expected_unstable_regions, regions_analyzed),
                max(expected_unstable_regions, regions_analyzed),
            ]

            deterministic_msi_variants = sum(
                len(
                    filter_variants_by_af_and_sample(
                        variants, af_threshold, sample_name
                    )
                )
                for variants in dp_ready_data["regions"].values()
            )
            msi_variants_impact = abs(
                expected_msi_variants - deterministic_msi_variants
            )
            msi_variants_range = [
                min(expected_msi_variants, deterministic_msi_variants),
                max(expected_msi_variants, deterministic_msi_variants),
            ]

            uncertain_regions = sample_uncertain_by_af.get(af_threshold, 0)

            sample_results[f"af_{af_threshold}"] = {
                # Core MSI Metrics
                "msi_score": round(msi_score, 2),
                "msi_score_deterministic": round(deterministic_msi_score, 2),
                "msi_score_range": [
                    round(msi_score_range[0], 2),
                    round(msi_score_range[1], 2),
                ],
                "analysis_impact": round(analysis_impact, 2),
                # Classification
                "msi_status": classify_msi_status(msi_score, msi_high_threshold),
                "msi_status_deterministic": classify_msi_status(
                    deterministic_msi_score, msi_high_threshold
                ),
                # Region Counts
                "total_regions": total_ms_regions,
                "unstable_regions": unstable_count,
                "deterministic_unstable_regions": regions_analyzed,
                "uncertain_regions": uncertain_regions,
                "probabilistic_stable_regions": total_ms_regions
                - unstable_count
                - uncertain_regions,
                "deterministic_stable_regions": total_ms_regions
                - regions_analyzed
                - uncertain_regions,
                # Expected Values with Uncertainty
                "expected_unstable_regions": round(expected_unstable_regions, 2),
                "expected_unstable_range": [
                    round(expected_unstable_range[0], 1),
                    round(expected_unstable_range[1], 1),
                ],
                "expected_unstable_impact": round(expected_unstable_impact, 1),
                "expected_variants_uncertainty": round(
                    expected_variants_uncertainty, 2
                ),
                "expected_unstable_uncertainty": round(
                    af_unstable_regions_uncertainty, 2
                ),
                "msi_score_statistical_uncertainty": round(
                    af_msi_score_statistical_uncertainty, 3
                ),
                # MSI Variants Analysis
                "expected_msi_variants": round(expected_msi_variants, 2),
                "deterministic_msi_variants": deterministic_msi_variants,
                "msi_variants_range": [
                    round(msi_variants_range[0], 1),
                    round(msi_variants_range[1], 1),
                ],
                "msi_variants_impact": round(msi_variants_impact, 1),
            }

        results[sample_name] = sample_results

        af_scores = [sample_results[f"af_{af}"]["msi_score"] for af in af_thresholds]
        print(
            f"[AF-EVOLUTION] {sample_name} evolution: {af_thresholds[0]}→{af_thresholds[-1]}: {af_scores[0]:.1f}%→{af_scores[-1]:.1f}%"
        )

    print(f"[AF-EVOLUTION] AF evolution analysis complete")
    return results


def test_af_evolution():
    """Test AF evolution with mock data including uncertain regions."""
    print("\n" + "=" * 60)
    print("TESTING AF EVOLUTION ANALYSIS WITH DEBUG")
    print("=" * 60)

    # Mock data with AF variations and uncertain regions
    mock_data = {
        "samples": ["tumor", "normal"],
        "af_thresholds": [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1, -2],
        "uncertain_by_af_by_sample": {
            "tumor": {
                1.0: 5,
                0.8: 8,
                0.6: 12,
                0.4: 15,
                0.2: 20,
                0.0: 25,
                -1: 30,
                -2: 2,
            },
            "normal": {
                1.0: 3,
                0.8: 6,
                0.6: 9,
                0.4: 12,
                0.2: 15,
                0.0: 20,
                -1: 25,
                -2: 1,
            },
        },
        "total_uncertain_regions": 35,
        "regions": {
            "chr1:1000-1020": [
                {
                    "dp_data": {
                        "p_present": 0.8,
                        "p_absent": 0.2,
                        "af_by_sample": {"tumor": 0.9, "normal": 0.1},
                    }
                },
                {
                    "dp_data": {
                        "p_present": 0.6,
                        "p_absent": 0.4,
                        "af_by_sample": {"tumor": 0.3, "normal": 0.05},
                    }
                },
            ],
            "chr1:2000-2020": [
                {
                    "dp_data": {
                        "p_present": 0.7,
                        "p_absent": 0.3,
                        "af_by_sample": {"tumor": 0.6, "normal": 0.02},
                    }
                }
            ],
            "chr1:3000-3030": [
                {
                    "dp_data": {
                        "p_present": 0.9,
                        "p_absent": 0.1,
                        "af_by_sample": {"tumor": 0.4, "normal": -1},
                    }
                }
            ],
            "chr1:4000-4030": [
                {
                    "dp_data": {
                        "p_present": 0.85,
                        "p_absent": 0.15,
                        "af_by_sample": {"tumor": -2, "normal": 0.8},
                    }
                }
            ],
        },
    }

    # Debug: Print mock data structure
    print("MOCK DATA DEBUG:")
    print(
        f"Uncertain regions for tumor at AF 0.8: {mock_data['uncertain_by_af_by_sample']['tumor'][0.8]}"
    )

    results = run_af_evolution_analysis(mock_data, total_ms_regions=1000)

    for sample_name, sample_data in results.items():
        print(f"\n{sample_name.upper()} AF EVOLUTION - DETAILED BREAKDOWN:")
        for af_key in sorted(sample_data.keys()):
            af_result = sample_data[af_key]
            print(f"\n  {af_key}:")
            print(
                f"    MSI Score: {af_result['msi_score']:.2f}% | Deterministic: {af_result['msi_score_deterministic']:.2f}%"
            )
            print(
                f"    Range: {af_result['msi_score_range'][0]:.1f}%-{af_result['msi_score_range'][1]:.1f}%"
            )
            print(
                f"    Analysis Impact: {af_result['analysis_impact']:.2f} percentage points"
            )
            print(
                f"    Classification: {af_result['msi_status']} (prob) vs {af_result['msi_status_deterministic']} (det)"
            )
            print(f"    Region Counts:")
            print(f"      Total: {af_result['total_regions']}")
            print(
                f"      Unstable (prob): {af_result['unstable_regions']} | Unstable (det): {af_result['deterministic_unstable_regions']}"
            )
            print(
                f"      Uncertain: {af_result['uncertain_regions']} | Stable (prob): {af_result['probabilistic_stable_regions']}"
            )
            print(f"    Expected Values:")
            print(
                f"      Unstable Regions: {af_result['expected_unstable_regions']:.1f} (Range: {af_result['expected_unstable_range'][0]:.1f}-{af_result['expected_unstable_range'][1]:.1f})"
            )
            print(
                f"      MSI Variants: {af_result['expected_msi_variants']:.1f} (Range: {af_result['msi_variants_range'][0]:.1f}-{af_result['msi_variants_range'][1]:.1f})"
            )
            print(
                f"      Variants Uncertainty: ± {af_result['expected_variants_uncertainty']:.2f}"
            )

    print("\n" + "=" * 60)
    print("AF EVOLUTION TESTING COMPLETE")
    print("=" * 60)


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


if __name__ == "__main__":
    test_af_evolution()
