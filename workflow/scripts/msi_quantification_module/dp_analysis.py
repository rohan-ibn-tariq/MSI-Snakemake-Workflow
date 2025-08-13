"""
MSI Quantification Module - Dynamic Programming Analysis
"""


def validate_probabilities(pp, pa, art):
    """
    Validate individual probability values upfront before processing.
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


def validate_af_value(af_value, sample_name):
    """Validate and normalize AF value for a specific sample."""
    if af_value is None:
        return -1, []
    
    try:
        af_float = float(af_value)
        if af_float < 0:
            return -1, [f"{sample_name}_negative_af_{af_float}"]
        if af_float > 1.0:
            return -1, [f"{sample_name}_af_exceeds_1.0_{af_float}"]
        return af_float, []
    except:
        return -1, [f"{sample_name}_invalid_af_type"]


def apply_4_step_imputation_to_variant(variant):
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
            "af_errors": []
        }
    }
    
    is_valid, prob_errors = validate_probabilities(pp, pa, art)
    if not is_valid:
        p_present = 0.5
        p_absent = 0.5
        audit_trail["imputation_method"] = "validation_fallback_uniform"
        audit_trail["calculation_steps"] = f"validation failed: {prob_errors}, applied uniform p_present=0.5, p_absent=0.5"
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
        
        if missing_count == 0:
            p_present = pp
            p_absent = pa + art
            audit_trail["imputation_method"] = "consolidation_direct"
            audit_trail["calculation_steps"] = f"p_absent = prob_absent + prob_artifact = {pa} + {art} = {p_absent}, p_present = {p_present}"
            
        elif missing_count == 1:
            if pp is None:
                derived = 1.0 - pa - art
                p_present = derived
                p_absent = pa + art
                audit_trail["imputation_method"] = "constraint_derived_present"
                audit_trail["calculation_steps"] = f"derived prob_present = 1.0 - {pa} - {art} = {derived}, p_absent = {pa} + {art} = {p_absent}"
                    
            elif pa is None:
                derived = 1.0 - pp - art
                p_present = pp
                p_absent = derived + art
                audit_trail["imputation_method"] = "constraint_derived_absent"
                audit_trail["calculation_steps"] = f"derived prob_absent = 1.0 - {pp} - {art} = {derived}, p_absent = {derived} + {art} = {p_absent}"
                    
            elif art is None:
                derived = 1.0 - pp - pa
                p_present = pp
                p_absent = pa + derived
                audit_trail["imputation_method"] = "constraint_derived_artifact"
                audit_trail["calculation_steps"] = f"derived prob_artifact = 1.0 - {pp} - {pa} = {derived}, p_absent = {pa} + {derived} = {p_absent}"
                    
        elif missing_count == 3:
            p_present = 0.5
            p_absent = 0.5
            audit_trail["imputation_method"] = "uniform_all_missing"
            audit_trail["calculation_steps"] = "all probabilities missing, applied uniform: present=0.5, absent=0.25, artifact=0.25, p_absent=0.5"
            
        else:
            known_value = pp if pp is not None else (pa if pa is not None else art)
            remaining = 1.0 - known_value
            
            if pp is None and pa is None:
                imputed_present = remaining * (2/3)
                imputed_absent = remaining * (1/3)
                p_present = imputed_present
                p_absent = imputed_absent + art
                audit_trail["imputation_method"] = "proportional_split_present_absent"
                audit_trail["calculation_steps"] = f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, absent = {remaining} * (1/3) = {imputed_absent}, p_absent = {imputed_absent} + {art} = {p_absent}"
                
            elif pp is None and art is None:
                imputed_present = remaining * (2/3)
                imputed_artifact = remaining * (1/3)
                p_present = imputed_present
                p_absent = pa + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_present_artifact"
                audit_trail["calculation_steps"] = f"remaining = {remaining}, present = {remaining} * (2/3) = {imputed_present}, artifact = {remaining} * (1/3) = {imputed_artifact}, p_absent = {pa} + {imputed_artifact} = {p_absent}"
                
            elif pa is None and art is None:
                imputed_absent = remaining * 0.5
                imputed_artifact = remaining * 0.5
                p_present = pp
                p_absent = imputed_absent + imputed_artifact
                audit_trail["imputation_method"] = "proportional_split_absent_artifact"
                audit_trail["calculation_steps"] = f"remaining = {remaining}, absent = {remaining} * 0.5 = {imputed_absent}, artifact = {remaining} * 0.5 = {imputed_artifact}, p_absent = {imputed_absent} + {imputed_artifact} = {p_absent}"
    
    sample_afs = variant.get("sample_afs", {})
    af_by_sample = {}
    
    for sample_name, af_value in sample_afs.items():
        normalized_af, af_errors = validate_af_value(af_value, sample_name)
        
        af_by_sample[sample_name] = normalized_af
        audit_trail["af_processing"]["af_errors"].extend(af_errors)
        
        if normalized_af == -1:
            audit_trail["af_processing"]["missing_af_samples"].append(sample_name)
            if af_value is None:
                audit_trail["af_processing"]["af_conversion"][sample_name] = "None -> -1"
            else:
                audit_trail["af_processing"]["af_conversion"][sample_name] = f"{af_value} -> -1 (error)"
        else:
            audit_trail["af_processing"]["af_conversion"][sample_name] = f"{af_value} -> {normalized_af}"
    
    if not af_by_sample:
        af_by_sample["sample"] = -1
        audit_trail["af_processing"]["missing_af_samples"].append("sample")
        audit_trail["af_processing"]["af_conversion"]["sample"] = "no_af_data -> -1"
    
    variant["dp_data"] = {
        "p_present": p_present,
        "p_absent": p_absent,
        "af_by_sample": af_by_sample,
        "audit_trail": audit_trail
    }


def prepare_variants_for_dp(results, imputation_method="uniform"):
    """Main entry point - complete Step 2 pipeline optimized into single function."""
    
    print(f"[DP-ENTRY] Starting Step 2 data preparation")
    print(f"[DP-ENTRY] Processing {len(results)} regions")
    
    # Extract variants and process in single pass
    region_variants = {}
    total_variants = 0
    sample_names_found = set()
    
    for region in results:
        region_id = region.get("region_id")
        variants = region.get("variants", [])
        
        if variants:            
            # Process each variant in this region
            for variant in variants:
                # Apply 4-step imputation
                apply_4_step_imputation_to_variant(variant)
                
                # Collect sample names
                dp_data = variant.get("dp_data", {})
                af_by_sample = dp_data.get("af_by_sample", {})
                sample_names_found.update(af_by_sample.keys())
                
                total_variants += 1
            
            # Store processed variants
            region_variants[region_id] = variants
    
    # Handle edge case
    if not region_variants:
        return {
            "error": "No variants found",
            "step": "data_preparation_failed",
            "total_regions": len(results)
        }
    
    # Create final structure
    sample_list = sorted(list(sample_names_found))
    af_thresholds = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1]
    
    dp_ready = {
        "af_thresholds": af_thresholds,
        "samples": sample_list,
        "regions": region_variants,
        "total_variants": total_variants,
        "total_regions": len(region_variants),
        "ready_for_dp": True,
        "entry_summary": {
            "input_regions": len(results),
            "regions_with_variants": len(region_variants),
            "total_variants": total_variants,
            "step_2_complete": True
        }
    }
    
    print(f"[DP-PREP] Processed {total_variants} variants in {len(region_variants)} regions")
    print(f"[DP-PREP] Found samples: {sample_list}")
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
            curr_col[i] = prev_col[i] * p_absent_k + prev_col[i-1] * p_present_k
        
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
    variance = sum((i - expected)**2 * prob for i, prob in enumerate(distribution))
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


def run_regional_msi_analysis(dp_ready_data, total_ms_regions_from_bed, unstable_threshold=0.5, msi_high_threshold=3.5):
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

    expected_unstable_regions = 0.0
    expected_msi_variants = 0.0  
    expected_unstable_uncertainties = []
    expected_variants_uncertainties = []

    
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
    
    msi_score = (unstable_count / total_ms_regions_from_bed) * 100 if total_ms_regions_from_bed > 0 else 0.0
    absolute_uncertainty, overall_uncertainty = propagate_uncertainties(region_uncertainties, total_ms_regions_from_bed)
    msi_status = classify_msi_status(msi_score, msi_high_threshold)
    
    expected_unstable_uncertainty = absolute_uncertainty
    expected_variants_uncertainty = absolute_uncertainty


    print(f"[REGIONAL-DP] Results: {unstable_count}/{total_ms_regions_from_bed} unstable regions")
    print(f"[REGIONAL-DP] MSI Score: {msi_score:.2f}% ± {overall_uncertainty:.2f}")
    print(f"[REGIONAL-DP] MSI Status: {msi_status}")
    print(f"[REGIONAL-DP] Total MS regions in BED: {total_ms_regions_from_bed}")
    print(f"[REGIONAL-DP] Regions with variants: {regions_with_variants}")
    print(f"[REGIONAL-DP] Regions without variants: {total_ms_regions_from_bed - regions_with_variants}")
    print(f"[REGIONAL-DP] Expected Unstable Regions: {expected_unstable_regions:.2f} ± {expected_unstable_uncertainty:.2f}")
    print(f"[REGIONAL-DP] Expected MSI Variants: {expected_msi_variants:.2f} ± {expected_variants_uncertainty:.2f}")

    return {
        "msi_score": round(msi_score, 2),
        "msi_uncertainty": round(overall_uncertainty, 2),
        "msi_status": msi_status,
        "unstable_regions": unstable_count,
        "total_regions": total_ms_regions_from_bed,
        "regions_with_variants": regions_with_variants,
        "regions_without_variants": total_ms_regions_from_bed - regions_with_variants,
        "expected_unstable_regions": round(expected_unstable_regions, 2),
        "expected_unstable_uncertainty": round(expected_unstable_uncertainty, 2),
        "expected_msi_variants": round(expected_msi_variants, 2),
        "expected_variants_uncertainty": round(expected_variants_uncertainty, 2),
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
        sample_af = af_by_sample.get(sample_name, -1)

        if af_threshold == -1:
            if sample_af == -1:
                filtered.append(variant)
        else:
            if sample_af >= af_threshold:
                filtered.append(variant)
    
    return filtered


def run_af_evolution_analysis(dp_ready_data, total_ms_regions, unstable_threshold=0.5):
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
        
        for af_threshold in af_thresholds:
            unstable_count = 0
            region_uncertainties = []
            regions_analyzed = 0
            
            for region_id, variants in dp_ready_data["regions"].items():
                filtered_variants = filter_variants_by_af_and_sample(variants, af_threshold, sample_name)
                
                if filtered_variants:
                    regions_analyzed += 1
                    
                    distribution = run_msi_dp(filtered_variants)
                    p_unstable = calculate_p_unstable(distribution)
                    uncertainty = calculate_std_dev(distribution)
                    
                    region_uncertainties.append(uncertainty)
                    
                    if p_unstable > unstable_threshold:
                        unstable_count += 1
            
            msi_score = (unstable_count / total_ms_regions) * 100 if total_ms_regions > 0 else 0.0
            msi_uncertainty = propagate_uncertainties(region_uncertainties, total_ms_regions)
            
            sample_results[f"af_{af_threshold}"] = {
                "msi_score": round(msi_score, 2),
                "msi_uncertainty": round(msi_uncertainty, 2),
                "unstable_regions": unstable_count,
                "regions_analyzed": regions_analyzed
            }
        
        results[sample_name] = sample_results
        
        af_scores = [sample_results[f"af_{af}"]["msi_score"] for af in af_thresholds]
        print(f"[AF-EVOLUTION] {sample_name} evolution: {af_thresholds[0]}→{af_thresholds[-1]}: {af_scores[0]:.1f}%→{af_scores[-1]:.1f}%")
    
    print(f"[AF-EVOLUTION] AF evolution analysis complete")
    return results


def test_af_evolution():
    """Test AF evolution with mock data."""
    print("\n" + "="*60)
    print("TESTING AF EVOLUTION ANALYSIS")
    print("="*60)
    
    # Mock data with AF variations
    mock_data = {
        "samples": ["tumor", "normal"],
        "af_thresholds": [1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -1],
        "regions": {
            "chr1:1000-1020": [
                {"dp_data": {"p_present": 0.8, "p_absent": 0.2, 
                            "af_by_sample": {"tumor": 0.9, "normal": 0.1}}},
                {"dp_data": {"p_present": 0.6, "p_absent": 0.4, 
                            "af_by_sample": {"tumor": 0.3, "normal": 0.05}}}
            ],
            "chr1:2000-2020": [
                {"dp_data": {"p_present": 0.7, "p_absent": 0.3, 
                            "af_by_sample": {"tumor": 0.6, "normal": 0.02}}}
            ],
            "chr1:3000-3030": [
                {"dp_data": {"p_present": 0.9, "p_absent": 0.1, 
                            "af_by_sample": {"tumor": 0.4, "normal": -1}}}  # Missing AF in normal
            ]
        }
    }
    
    results = run_af_evolution_analysis(mock_data, total_ms_regions=1000)
    
    for sample_name, sample_data in results.items():
        print(f"\n{sample_name.upper()} AF EVOLUTION:")
        for af_key in sorted(sample_data.keys()):
            af_result = sample_data[af_key]
            print(f"  {af_key}: {af_result['msi_score']:.2f}% ± {af_result['msi_uncertainty']:.2f}% "
                  f"({af_result['unstable_regions']} unstable, {af_result['regions_analyzed']} analyzed)")
    
    print("\n" + "="*60)
    print("AF EVOLUTION TESTING COMPLETE")
    print("="*60)


if __name__ == "__main__":
    test_af_evolution()
