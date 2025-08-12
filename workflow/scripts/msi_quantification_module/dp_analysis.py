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
            print(f"[DP-PREP] Processing region {region_id} with {len(variants)} variants")
            
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
    print("[DP-ENTRY] Data Preparation for DP complete")

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

