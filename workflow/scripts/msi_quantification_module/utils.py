"""
MSI Quantification Module - Utility Functions

Common utility functions for MSI analysis including probability conversions
and data processing helpers.
"""

import math


def phred_to_prob(phred_score):
    """
    Convert PHRED score to probability following VarLociRaptor conventions.
    
    Args:
        phred_score: PHRED quality score (float, None, or inf)
        
    Returns:
        float or None: Probability value
        - None: Missing data
        - 0.0: Perfect confidence (infinite PHRED)  
        - 0.0-1.0: Converted probability for normal PHRED scores
    """
    if phred_score is None:
        return None  # Missing data - preserve for analysis-level decisions
    elif math.isinf(phred_score):
        return 0.0   # Perfect confidence = 0.0 error probability
    else:
        return 10 ** (-float(phred_score) / 10)
