"""Cross-validation of BMD/BCA results between the Hologic APEX XML tag
(0019,1000) and the companion DICOM Structured Report.

The two sources use different region-name conventions:
* XML  — short English labels: "L Arm", "R Arm", "Trunk", "Total", …
* SR   — anatomic codes whose CodeMeaning varies by locale/APEX version.

This module normalises both to a shared canonical vocabulary and produces a
tidy comparison DataFrame with absolute and relative deltas for every matched
quantity.

Typical usage
-------------
    from hologic_dxa.dicom.bmd_comparison import compare_apex_xml_vs_sr
    import pydicom
    ds_scan = pydicom.dcmread("scan.dcm", force=True)
    ds_sr   = pydicom.dcmread("sr.dcm",   force=True)

    from hologic_dxa.dicom.hologic_xml      import parse_apex_xml
    from hologic_dxa.dicom.structured_report import parse_structured_report
    df = compare_apex_xml_vs_sr(parse_apex_xml(ds_scan), parse_structured_report(ds_sr))
    print(df.to_string())
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

from hologic_dxa.dicom.hologic_xml import ApexResults
from hologic_dxa.dicom.structured_report import SRMeasurement


# ---------------------------------------------------------------------------
# Region-name normalisation
# ---------------------------------------------------------------------------

# Canonical region names (XML vocabulary is the reference).
# Each entry: (canonical, [aliases…])
_REGION_ALIASES: list[tuple[str, list[str]]] = [
    ("L Arm",     ["left arm", "l arm", "arm left", "left upper extremity"]),
    ("R Arm",     ["right arm", "r arm", "arm right", "right upper extremity"]),
    ("L Ribs",   ["left ribs", "l ribs", "ribs left", "left rib"]),
    ("R Ribs",   ["right ribs", "r ribs", "ribs right", "right rib"]),
    ("T Spine",  ["thoracic spine", "t spine", "thoracic", "t-spine"]),
    ("L Spine",  ["lumbar spine", "l spine", "lumbar", "l-spine", "l1-l4", "l1l4"]),
    ("Pelvis",   ["pelvis", "pelvic", "hip"]),
    ("Trunk",    ["trunk", "torso"]),
    ("L Leg",    ["left leg", "l leg", "leg left", "left lower extremity"]),
    ("R Leg",    ["right leg", "r leg", "leg right", "right lower extremity"]),
    ("Subtotal", ["subtotal", "sub total", "sub-total"]),
    ("Head",     ["head", "skull", "cranium"]),
    ("Total",    ["total", "whole body", "total body"]),
]

# Build lookup: normalised alias → canonical
_ALIAS_TO_CANONICAL: dict[str, str] = {}
for _canon, _aliases in _REGION_ALIASES:
    _ALIAS_TO_CANONICAL[_canon.lower()] = _canon
    for _a in _aliases:
        _ALIAS_TO_CANONICAL[_a.lower()] = _canon


def normalise_region(name: str) -> str:
    """Return the canonical region name, or the stripped input if unknown."""
    key = re.sub(r"\s+", " ", name.strip()).lower()
    return _ALIAS_TO_CANONICAL.get(key, name.strip())


# ---------------------------------------------------------------------------
# SR quantity mapping
# ---------------------------------------------------------------------------

# DICOM concept CodeValues that correspond to XML columns.
# Hologic APEX 13.x uses these DCM + 99SDM codes; extend as needed.
_SR_CODE_TO_COLUMN: dict[str, str] = {
    # DCM codes (confirmed in Hologic SR)
    "122490": "bmc_g",          # Bone Mineral Content
    "122492": "bmd_g_cm2",      # Bone Mineral Density
    "122488": "area_cm2",       # Area
    # Proprietary / extended codes (may vary; kept as fallback)
    "122489": "lean_g",         # Lean Tissue Mass
    "122491": "fat_g",          # Fat Tissue Mass
    "122493": "fat_percent",    # Percent Fat
    "122495": "total_g",        # Total Tissue Mass
    "122496": "lean_bmc_g",     # Lean + BMC
}

# Code-meaning fallback (when code_value is missing / non-standard)
_SR_MEANING_TO_COLUMN: dict[str, str] = {
    "bone mineral content":   "bmc_g",
    "bone mineral density":   "bmd_g_cm2",
    "area":                   "area_cm2",
    "lean tissue mass":       "lean_g",
    "fat tissue mass":        "fat_g",
    "percent fat":            "fat_percent",
    "total tissue mass":      "total_g",
    "lean+bmc":               "lean_bmc_g",
    "lean + bmc":             "lean_bmc_g",
}


def _sr_to_column(m: SRMeasurement) -> str | None:
    """Map an SRMeasurement to its XML column name, or None if unmapped."""
    col = _SR_CODE_TO_COLUMN.get(m.code_value or "")
    if col:
        return col
    meaning = (m.code_meaning or "").strip().lower()
    return _SR_MEANING_TO_COLUMN.get(meaning)


# ---------------------------------------------------------------------------
# Core comparison
# ---------------------------------------------------------------------------

#: Columns that exist in both XML and (potentially) SR.
COMPARABLE_COLUMNS = [
    "area_cm2", "bmc_g", "bmd_g_cm2",
    "fat_g", "lean_g", "lean_bmc_g", "total_g", "fat_percent",
]


@dataclass
class ComparisonResult:
    """Outcome of cross-validating XML vs SR for one study."""

    df: pd.DataFrame
    """Tidy DataFrame — one row per (region, quantity)."""

    unmatched_xml_regions: list[str] = field(default_factory=list)
    """XML regions that had no SR counterpart."""

    unmatched_sr_regions: list[str] = field(default_factory=list)
    """SR regions that had no XML counterpart."""

    unmapped_sr_codes: list[tuple[str, str]] = field(default_factory=list)
    """(code_value, code_meaning) pairs from SR that were not mapped to a column."""


def compare_apex_xml_vs_sr(
    xml_results: ApexResults,
    sr_measurements: list[SRMeasurement],
) -> ComparisonResult:
    """Cross-validate APEX XML vs SR measurements.

    Parameters
    ----------
    xml_results:
        Parsed output of :func:`~hologic_dxa.dicom.hologic_xml.parse_apex_xml`.
    sr_measurements:
        Parsed output of :func:`~hologic_dxa.dicom.structured_report.parse_structured_report`.

    Returns
    -------
    ComparisonResult
        ``.df`` has columns: region, quantity, xml_value, sr_value,
        abs_diff (xml − sr), rel_diff_pct (100 × abs_diff / sr_value).
    """
    # ── 1. Build XML wide table ────────────────────────────────────────────
    xml_rows: list[dict] = []
    for m in xml_results.measurements:
        region = normalise_region(m.region)
        for col in COMPARABLE_COLUMNS:
            val = getattr(m, col, None)
            if val is not None:
                xml_rows.append({"region": region, "quantity": col, "xml_value": val})

    xml_df = pd.DataFrame(xml_rows, columns=["region", "quantity", "xml_value"])

    # ── 2. Build SR wide table ─────────────────────────────────────────────
    unmapped: list[tuple[str, str]] = []
    sr_rows: list[dict] = []
    for m in sr_measurements:
        col = _sr_to_column(m)
        if col is None:
            unmapped.append((m.code_value or "", m.code_meaning or ""))
            continue
        if m.value is None:
            continue
        region_raw = m.anatomic_region_meaning or m.tracking_id or ""
        region = normalise_region(region_raw)
        sr_rows.append({"region": region, "quantity": col, "sr_value": m.value})

    sr_df = pd.DataFrame(sr_rows, columns=["region", "quantity", "sr_value"])

    # ── 3. Join ────────────────────────────────────────────────────────────
    merged = pd.merge(xml_df, sr_df, on=["region", "quantity"], how="outer")

    # ── 4. Deltas ──────────────────────────────────────────────────────────
    merged["abs_diff"] = merged["xml_value"] - merged["sr_value"]
    merged["rel_diff_pct"] = (
        100.0 * merged["abs_diff"] / merged["sr_value"].replace(0, float("nan"))
    )
    merged = merged.sort_values(["region", "quantity"]).reset_index(drop=True)

    # ── 5. Unmatched sets ─────────────────────────────────────────────────
    xml_regions = set(xml_df["region"].unique())
    sr_regions  = set(sr_df["region"].unique())
    unmatched_xml = sorted(xml_regions - sr_regions)
    unmatched_sr  = sorted(sr_regions  - xml_regions)

    return ComparisonResult(
        df=merged,
        unmatched_xml_regions=unmatched_xml,
        unmatched_sr_regions=unmatched_sr,
        unmapped_sr_codes=sorted(set(unmapped)),
    )


def comparison_summary(result: ComparisonResult) -> pd.DataFrame:
    """Summarise per-quantity agreement statistics across all matched regions.

    Returns a DataFrame indexed by quantity with columns:
    n_pairs, mean_abs_diff, max_abs_diff, mean_rel_diff_pct, max_rel_diff_pct.
    """
    matched = result.df.dropna(subset=["xml_value", "sr_value"])
    if matched.empty:
        return pd.DataFrame(
            columns=["quantity", "n_pairs", "mean_abs_diff", "max_abs_diff",
                     "mean_rel_diff_pct", "max_rel_diff_pct"]
        )
    grp = matched.groupby("quantity")
    summary = pd.DataFrame({
        "n_pairs":          grp["abs_diff"].count(),
        "mean_abs_diff":    grp["abs_diff"].apply(lambda x: x.abs().mean()),
        "max_abs_diff":     grp["abs_diff"].apply(lambda x: x.abs().max()),
        "mean_rel_diff_pct": grp["rel_diff_pct"].apply(lambda x: x.abs().mean()),
        "max_rel_diff_pct":  grp["rel_diff_pct"].apply(lambda x: x.abs().max()),
    }).reset_index()
    return summary
