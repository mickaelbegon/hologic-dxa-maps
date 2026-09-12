"""Regional comparison between integrated map masses and APEX SR values."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class RegionalMassComparison:
    """Comparison of integrated vs. reference mass for one region and component."""

    region: str
    component: str        # "fat", "lean", "bmc", "total"
    integrated_g: float
    reference_g: float
    absolute_error_g: float
    relative_error: float
    bland_altman_mean_g: float
    passed: bool


def compare_regional_masses(
    integrated: dict[str, dict[str, float]],
    sr_df: pd.DataFrame,
    tolerance_relative: float = 0.05,
    tolerance_absolute_g: float = 50.0,
) -> list[RegionalMassComparison]:
    """Compare integrated masses against SR reference values.

    Parameters
    ----------
    integrated:
        {region: {component: mass_g}} from compute_mass_properties()
    sr_df:
        Long-format SR DataFrame with columns: canonical_name, value, unit
    tolerance_relative, tolerance_absolute_g:
        PROVISIONAL thresholds. A comparison passes if EITHER is met.

    Returns
    -------
    List of RegionalMassComparison, one per (region, component) pair found.
    """
    results = []

    for region, props in integrated.items():
        for component in ("fat", "lean", "bmc", "total"):
            key = f"{component}_g"
            if key not in props:
                continue
            integrated_g = props[key]

            # Look up reference from SR
            ref_val = _find_sr_reference(sr_df, region, component)
            if ref_val is None:
                continue

            abs_err = abs(integrated_g - ref_val)
            rel_err = (integrated_g - ref_val) / ref_val if ref_val != 0 else float("inf")
            ba_mean = (integrated_g + ref_val) / 2.0

            passed = (
                abs(rel_err) <= tolerance_relative
                or abs_err <= tolerance_absolute_g
            )

            results.append(RegionalMassComparison(
                region=region,
                component=component,
                integrated_g=integrated_g,
                reference_g=ref_val,
                absolute_error_g=abs_err,
                relative_error=rel_err,
                bland_altman_mean_g=ba_mean,
                passed=passed,
            ))

    return results


def _find_sr_reference(
    sr_df: pd.DataFrame,
    region: str,
    component: str,
) -> float | None:
    """Locate a reference mass value in the SR DataFrame.

    This uses the canonical_name column from the concept mapping normalization.
    If the SR has not been normalized, this function returns None.
    """
    if sr_df is None or sr_df.empty:
        return None
    if "canonical_name" not in sr_df.columns:
        return None

    # Map component to canonical name fragment
    _COMPONENT_MAP = {
        "fat": "fat_mass",
        "lean": "lean_mass",
        "bmc": "bone_mineral_content",
        "total": "total_mass",
    }
    target = _COMPONENT_MAP.get(component)
    if target is None:
        return None

    region_lower = region.lower().replace(" ", "_")
    mask = (
        sr_df["canonical_name"].str.contains(target, case=False, na=False)
        & sr_df["canonical_name"].str.contains(region_lower, case=False, na=False)
    )
    matches = sr_df[mask]
    if matches.empty:
        return None

    val = matches["value"].iloc[0]
    unit = str(matches.get("unit_meaning", pd.Series(["g"])).iloc[0]).lower()

    # Convert to grams if needed
    if "kg" in unit:
        return float(val) * 1000.0
    if "mg" in unit:
        return float(val) / 1000.0
    return float(val)
