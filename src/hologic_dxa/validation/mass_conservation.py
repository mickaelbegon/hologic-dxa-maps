"""Mass conservation validation for quantitative DXA pixel maps.

Compares pixel-integrated regional masses against reference values obtained
from Hologic APEX Structured Reports (SR) or APEX CSV exports.

IMPORTANT — threshold provenance
---------------------------------
All pass/fail thresholds used here are PROVISIONAL.  They have NOT been
validated against phantom data and must be calibrated empirically before any
result can be considered a scientific validation.  Any ``passed=True`` result
produced by this module is conditional on proper calibration.

Expected SR DataFrame format (long-form)
-----------------------------------------
Columns: ``region``, ``component``, ``value_g``

- ``region``: region name string (e.g. ``"whole_body"``, ``"left_arm"``).
- ``component``: one of ``"fat_g"``, ``"lean_g"``, ``"bmc_g"``, ``"total_g"``.
- ``value_g``: reference mass in grams (float).

Expected APEX CSV format
-------------------------
Columns: ``region``, ``fat_g``, ``lean_g``, ``bmc_g``, ``total_g``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from hologic_dxa import __version__
from hologic_dxa.logging import get_logger
from hologic_dxa.maps.integration import compute_mass_properties

if TYPE_CHECKING:
    from hologic_dxa.models import QuantitativeMapBundle

log = get_logger(__name__)


@dataclass
class RegionComparison:
    """Comparison of integrated vs reference mass for one region/component pair."""

    region_name: str
    integrated_g: float
    reference_g: float
    absolute_error_g: float
    relative_error: float
    n_valid_pixels: int
    passed: bool


@dataclass
class ValidationReport:
    """Full mass conservation validation report."""

    timestamp_utc: datetime
    pipeline_version: str
    comparisons: list[RegionComparison]
    overall_passed: bool
    thresholds_used: dict
    warnings: list[str]
    provenance: dict

    def to_dataframe(self) -> pd.DataFrame:
        """Return comparisons as a pandas DataFrame."""
        return pd.DataFrame(
            [
                {
                    "region_name": c.region_name,
                    "integrated_g": c.integrated_g,
                    "reference_g": c.reference_g,
                    "absolute_error_g": c.absolute_error_g,
                    "relative_error": c.relative_error,
                    "n_valid_pixels": c.n_valid_pixels,
                    "passed": c.passed,
                }
                for c in self.comparisons
            ]
        )

    def to_json(self) -> str:
        """Serialize the full report to a JSON string."""
        d = {
            "timestamp_utc": self.timestamp_utc.isoformat(),
            "pipeline_version": self.pipeline_version,
            "overall_passed": self.overall_passed,
            "thresholds_used": self.thresholds_used,
            "warnings": self.warnings,
            "provenance": self.provenance,
            "comparisons": [
                {
                    "region_name": c.region_name,
                    "integrated_g": c.integrated_g,
                    "reference_g": c.reference_g,
                    "absolute_error_g": c.absolute_error_g,
                    "relative_error": c.relative_error,
                    "n_valid_pixels": c.n_valid_pixels,
                    "passed": c.passed,
                }
                for c in self.comparisons
            ],
        }
        return json.dumps(d, indent=2)


def validate_mass_conservation(
    bundle: "QuantitativeMapBundle",
    sr_results: pd.DataFrame | None = None,
    apex_csv: pd.DataFrame | None = None,
    mass_tolerance_relative: float = 0.05,
    mass_tolerance_absolute_g: float = 50.0,
) -> ValidationReport:
    """Compare integrated masses from pixel maps vs SR/APEX reference values.

    Parameters
    ----------
    bundle:
        ``QuantitativeMapBundle`` with maps, masks, and geometry.
    sr_results:
        Long-format SR DataFrame with columns ``region``, ``component``,
        ``value_g``.  Expected components: ``fat_g``, ``lean_g``, ``bmc_g``,
        ``total_g``.
    apex_csv:
        Wide-format APEX DataFrame with columns ``region``, ``fat_g``,
        ``lean_g``, ``bmc_g``, ``total_g``.
    mass_tolerance_relative:
        PROVISIONAL maximum acceptable |relative error| = |integrated - reference| /
        reference.  Must be established empirically using phantom data.
    mass_tolerance_absolute_g:
        PROVISIONAL maximum acceptable absolute error in grams.  Must be
        established empirically using phantom data.

    Returns
    -------
    ValidationReport
        Populated report.  When no reference data is supplied all comparisons
        are skipped and a warning is added.

    Notes
    -----
    - Thresholds are PROVISIONAL and must be established empirically using
      phantom data.
    - Any ``passed=True`` result is conditional on calibration quality.
    - If no reference data is provided, the report contains no comparisons
      and a warning stating that validation was skipped.
    """
    timestamp = datetime.now(tz=timezone.utc)
    thresholds: dict = {
        "mass_tolerance_relative": mass_tolerance_relative,
        "mass_tolerance_absolute_g": mass_tolerance_absolute_g,
        "note": (
            "PROVISIONAL — not validated against phantom data. "
            "Calibrate empirically before treating pass/fail as scientifically meaningful."
        ),
    }
    report_warnings: list[str] = []

    if sr_results is None and apex_csv is None:
        log.warning(
            "validation_skipped",
            reason="no_reference_data",
        )
        report_warnings.append(
            "No reference data supplied (sr_results and apex_csv are both None). "
            "Mass conservation validation was skipped. "
            "Provide SR or APEX CSV data to enable comparison."
        )
        return ValidationReport(
            timestamp_utc=timestamp,
            pipeline_version=__version__,
            comparisons=[],
            overall_passed=False,
            thresholds_used=thresholds,
            warnings=report_warnings,
            provenance=_extract_provenance(bundle),
        )

    # Compute integrated masses from pixel maps
    try:
        mass_props = compute_mass_properties(bundle)
    except ValueError as exc:
        raise ValueError(
            f"Cannot compute integrated masses from bundle: {exc}"
        ) from exc

    # Build reference lookup: {(region, component): value_g}
    reference: dict[tuple[str, str], float] = {}

    if sr_results is not None:
        reference.update(_parse_sr_long(sr_results, report_warnings))

    if apex_csv is not None:
        reference.update(_parse_apex_wide(apex_csv, report_warnings))

    comparisons: list[RegionComparison] = []
    components = ("fat_g", "lean_g", "bmc_g", "total_g")

    for region_name, props in mass_props.items():
        valid_mask = _get_valid_mask(bundle, region_name)
        n_valid = int(np.sum(valid_mask)) if valid_mask is not None else -1

        for component in components:
            integrated = props.get(component, float("nan"))
            if np.isnan(integrated):
                continue

            key = (region_name, component)
            if key not in reference:
                log.debug(
                    "reference_not_found",
                    region=region_name,
                    component=component,
                )
                continue

            ref_val = reference[key]
            abs_err = integrated - ref_val
            rel_err = abs_err / ref_val if ref_val != 0.0 else float("nan")
            passed = (
                abs(abs_err) <= mass_tolerance_absolute_g
                and (np.isnan(rel_err) or abs(rel_err) <= mass_tolerance_relative)
            )

            comparisons.append(
                RegionComparison(
                    region_name=f"{region_name}/{component}",
                    integrated_g=integrated,
                    reference_g=ref_val,
                    absolute_error_g=abs_err,
                    relative_error=rel_err,
                    n_valid_pixels=n_valid,
                    passed=passed,
                )
            )

    if not comparisons:
        report_warnings.append(
            "No region/component pairs could be matched between integrated maps and "
            "reference data.  Check that region names and component names align."
        )

    overall_passed = bool(comparisons) and all(c.passed for c in comparisons)

    log.info(
        "mass_conservation_validated",
        n_comparisons=len(comparisons),
        overall_passed=overall_passed,
    )

    return ValidationReport(
        timestamp_utc=timestamp,
        pipeline_version=__version__,
        comparisons=comparisons,
        overall_passed=overall_passed,
        thresholds_used=thresholds,
        warnings=report_warnings,
        provenance=_extract_provenance(bundle),
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_sr_long(
    df: pd.DataFrame,
    warnings: list[str],
) -> dict[tuple[str, str], float]:
    """Parse a long-format SR DataFrame into a {(region, component): value_g} dict."""
    required = {"region", "component", "value_g"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        warnings.append(
            f"SR DataFrame is missing expected columns: {sorted(missing_cols)}. "
            "Expected columns are 'region', 'component', 'value_g'. "
            "SR results were not used in the comparison."
        )
        return {}

    result: dict[tuple[str, str], float] = {}
    for _, row in df.iterrows():
        region = str(row["region"])
        component = str(row["component"])
        try:
            value = float(row["value_g"])
        except (TypeError, ValueError):
            warnings.append(
                f"SR row (region='{region}', component='{component}') has a "
                f"non-numeric value_g ('{row['value_g']}'); skipped."
            )
            continue
        result[(region, component)] = value
    return result


def _parse_apex_wide(
    df: pd.DataFrame,
    warnings: list[str],
) -> dict[tuple[str, str], float]:
    """Parse a wide-format APEX DataFrame into a {(region, component): value_g} dict."""
    if "region" not in df.columns:
        warnings.append(
            "APEX CSV DataFrame is missing 'region' column. "
            "APEX data were not used in the comparison."
        )
        return {}

    component_cols = [c for c in ("fat_g", "lean_g", "bmc_g", "total_g") if c in df.columns]
    if not component_cols:
        warnings.append(
            "APEX CSV DataFrame contains none of the expected component columns "
            "('fat_g', 'lean_g', 'bmc_g', 'total_g'). "
            "APEX data were not used in the comparison."
        )
        return {}

    result: dict[tuple[str, str], float] = {}
    for _, row in df.iterrows():
        region = str(row["region"])
        for comp in component_cols:
            try:
                value = float(row[comp])
            except (TypeError, ValueError):
                warnings.append(
                    f"APEX row (region='{region}', component='{comp}') has a "
                    f"non-numeric value ('{row[comp]}'); skipped."
                )
                continue
            result[(region, comp)] = value
    return result


def _extract_provenance(bundle: "QuantitativeMapBundle") -> dict:
    """Extract provenance info from bundle for the report."""
    prov = getattr(bundle, "provenance", None)
    if prov is None:
        return {}
    if hasattr(prov, "to_dict"):
        return prov.to_dict()
    if isinstance(prov, dict):
        return prov
    return {}


def _get_valid_mask(
    bundle: "QuantitativeMapBundle", region_name: str
) -> np.ndarray | None:
    """Return the valid mask for a region, if available."""
    masks = getattr(bundle, "masks", None)
    if masks is None:
        return None
    if region_name == "whole_body":
        return getattr(masks, "valid", None) or getattr(masks, "body", None)
    return None
