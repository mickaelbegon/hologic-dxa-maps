"""Repeatability and reproducibility analysis for DXA quantitative maps.

Quantifies measurement variability from:
- Test-retest acquisitions (same subject, repositioned)
- Between-session acquisitions
- Cross-APEX-version comparisons

All metrics follow ISO 5725 and ISCD guidelines for DXA precision assessment.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def coefficient_of_variation(values: list[float]) -> float:
    """CV% = (std / mean) × 100. Requires ≥ 2 values."""
    arr = np.array(values, dtype=np.float64)
    if len(arr) < 2:
        raise ValueError("At least 2 values are required to compute CV.")
    return float(np.std(arr, ddof=1) / np.mean(arr) * 100)


def repeatability_coefficient(values: list[float]) -> float:
    """RC = 2.77 × SD (95% limit of agreement for paired measurements, ISO 5725)."""
    arr = np.array(values, dtype=np.float64)
    if len(arr) < 2:
        raise ValueError("At least 2 values are required.")
    return float(2.77 * np.std(arr, ddof=1))


def intraclass_correlation(
    measurements: pd.DataFrame,
    subject_col: str = "research_id",
    value_col: str = "mass_g",
    session_col: str = "session",
) -> dict[str, float]:
    """Two-way mixed ICC(2,1) for test-retest reliability.

    Requires a DataFrame with one row per (subject, session) pair.
    Returns ICC estimate, 95% CI, and within-subject CV.

    Raises
    ------
    NotImplementedError
        Until ICC computation is validated with real test-retest data.
    """
    raise NotImplementedError(
        "ICC computation is not yet implemented. "
        "Requires at least 10 subject pairs with 2 sessions each "
        "to provide stable estimates. "
        "Implement using the pingouin library (pip install pingouin) "
        "when test-retest data are available."
    )
