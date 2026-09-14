"""Hologic APEX XML results parser — private tag (0019,1000).

Hologic APEX embeds a proprietary XML + JavaScript document in the private
DICOM tag (0019,1000).  This module extracts region-level body-composition
and bone-density measurements from that document without relying on any
Hologic SDK or on the DICOM Structured Report.

Format notes
------------
* The tag contains UTF-8-encoded text mixing XML markup with embedded
  ``<JavaScriptFile>`` sections that hold variable assignments.
* The results table is encoded as::

      ResultsTable1[ row][ col] = "value";

  where ``row=0`` is the header row (column labels, possibly with HTML tags)
  and subsequent rows hold per-region values.
* Numeric values may be empty strings for regions where a quantity is not
  computed (e.g. BMD is not reported for the Trunk as a whole).
* The ``TableNote1`` assignment carries global calibration metadata including
  the ACF (Air Calibration Factor) and BCF (Bone Calibration Factor).

Security constraints
--------------------
* The raw XML text is treated as untrusted data.  Only numeric and short
  string fields are extracted; no code from the embedded JavaScript is
  evaluated.
* Returns empty/None results rather than raising when the tag is absent or
  malformed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd
import pydicom
import structlog

log = structlog.get_logger(__name__)

# Private tag (0019,1000) — Hologic private creator "HOLOGIC" in group 0019
_HOLOGIC_XML_TAG = (0x0019, 0x1000)

# Regex: ResultsTable1[ row][ col] = "value";   (spaces inside brackets ok)
_TABLE_RE = re.compile(
    r'ResultsTable1\[\s*(\d+)\s*\]\[\s*(\d+)\s*\]\s*=\s*"([^"]*)"'
)
# Regex: strip HTML tags from values/labels
_HTML_TAG_RE = re.compile(r"<[^>]+>")

# Regex for global scalar assignments in JavaScript blocks, e.g.  Height = "165.0 cm";
_JS_VAR_RE = re.compile(r'^(\w+)\s*=\s*"([^"]*)"\s*;', re.MULTILINE)

# Pattern to pull ACF / BCF out of the TableNote string
_ACF_RE = re.compile(r"ACF\s*=\s*([\d.]+)")
_BCF_RE = re.compile(r"BCF\s*=\s*([\d.]+)")


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class ApexRegionMeasurement:
    """Body-composition and bone-density measurements for one body region."""

    region: str
    area_cm2: float | None = None
    bmc_g: float | None = None
    bmd_g_cm2: float | None = None
    fat_g: float | None = None
    lean_g: float | None = None
    lean_bmc_g: float | None = None
    total_g: float | None = None
    fat_percent: float | None = None
    tscore: float | None = None
    zscore: float | None = None


@dataclass
class ApexScanInfo:
    """Scan and patient metadata from the Hologic XML header."""

    height_cm: float | None = None
    weight_kg: float | None = None
    age: int | None = None
    sex: str | None = None
    scan_date: str | None = None
    scan_protocol: str | None = None
    scan_mode: str | None = None
    acf: float | None = None
    bcf: float | None = None
    calibration_note: str | None = None


@dataclass
class ApexResults:
    """Complete results parsed from the Hologic APEX XML tag (0019,1000)."""

    measurements: list[ApexRegionMeasurement] = field(default_factory=list)
    scan_info: ApexScanInfo = field(default_factory=ApexScanInfo)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_apex_xml(ds: pydicom.Dataset) -> ApexResults:
    """Extract APEX body-composition results from a Hologic DICOM dataset.

    Parameters
    ----------
    ds:
        A pydicom Dataset.  May or may not contain the Hologic XML tag.

    Returns
    -------
    ApexResults
        Parsed results.  ``measurements`` is empty when the tag is absent or
        contains no recognisable table entries.  Anomalies are recorded in
        ``warnings`` rather than raised as exceptions.
    """
    result = ApexResults()

    if _HOLOGIC_XML_TAG not in ds:
        log.debug("apex_xml_tag_absent", tag=hex(_HOLOGIC_XML_TAG[0]))
        return result

    raw_bytes = ds[_HOLOGIC_XML_TAG].value
    try:
        xml_text = (
            bytes(raw_bytes).decode("utf-8", errors="replace")
            if isinstance(raw_bytes, (bytes, bytearray, memoryview))
            else str(raw_bytes)
        )
    except Exception as exc:  # noqa: BLE001
        result.warnings.append(f"Could not decode (0019,1000) as UTF-8: {exc}")
        return result

    _parse_results_table(xml_text, result)
    _parse_scan_info(xml_text, result)
    log.info(
        "apex_xml_parsed",
        region_count=len(result.measurements),
        warning_count=len(result.warnings),
    )
    return result


def apex_results_to_dataframe(results: ApexResults) -> pd.DataFrame:
    """Convert :class:`ApexResults` to a tidy long-format DataFrame.

    Each row is one body region.  Returns an empty DataFrame (correct columns)
    when ``results.measurements`` is empty.
    """
    _cols = [
        "region", "area_cm2", "bmc_g", "bmd_g_cm2",
        "fat_g", "lean_g", "lean_bmc_g", "total_g",
        "fat_percent", "tscore", "zscore",
    ]
    if not results.measurements:
        return pd.DataFrame(columns=_cols)

    rows = [
        {
            "region": m.region,
            "area_cm2": m.area_cm2,
            "bmc_g": m.bmc_g,
            "bmd_g_cm2": m.bmd_g_cm2,
            "fat_g": m.fat_g,
            "lean_g": m.lean_g,
            "lean_bmc_g": m.lean_bmc_g,
            "total_g": m.total_g,
            "fat_percent": m.fat_percent,
            "tscore": m.tscore,
            "zscore": m.zscore,
        }
        for m in results.measurements
    ]
    return pd.DataFrame(rows, columns=_cols)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _strip_html(text: str) -> str:
    return _HTML_TAG_RE.sub("", text).strip()


def _to_float(raw: str) -> float | None:
    """Parse a numeric string, returning None for empty or unparseable input."""
    s = _strip_html(raw).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_results_table(xml_text: str, result: ApexResults) -> None:
    """Fill ``result.measurements`` from ResultsTable1 entries."""
    # Collect all (row, col, value) entries
    raw_table: dict[int, dict[int, str]] = {}
    for match in _TABLE_RE.finditer(xml_text):
        row, col, val = int(match.group(1)), int(match.group(2)), match.group(3)
        raw_table.setdefault(row, {})[col] = val

    if not raw_table:
        log.debug("apex_xml_no_table_entries")
        return

    # Row 0 = header; rows 1+ = data
    n_rows = max(raw_table.keys())
    if n_rows < 1:
        result.warnings.append("ResultsTable1 has no data rows (only header).")
        return

    # Column index → field mapping (0-indexed, col 0 = region label)
    # Layout confirmed for APEX 13.6; earlier versions may differ.
    _COL_MAP = {
        1: "area_cm2",
        2: "bmc_g",
        3: "bmd_g_cm2",
        4: "fat_g",
        5: "lean_g",
        6: "lean_bmc_g",
        7: "total_g",
        8: "fat_percent",
        9: "tscore",
        10: "zscore",
    }

    for row_idx in sorted(raw_table.keys()):
        if row_idx == 0:
            continue  # skip header
        row_data = raw_table[row_idx]
        region_raw = row_data.get(0, "")
        region = _strip_html(region_raw)
        if not region:
            result.warnings.append(f"Row {row_idx} has empty region label — skipped.")
            continue

        m = ApexRegionMeasurement(region=region)
        for col_idx, attr_name in _COL_MAP.items():
            raw_val = row_data.get(col_idx, "")
            setattr(m, attr_name, _to_float(raw_val))

        result.measurements.append(m)


def _parse_scan_info(xml_text: str, result: ApexResults) -> None:
    """Fill ``result.scan_info`` from JavaScript variable assignments."""
    info = result.scan_info
    js_vars: dict[str, str] = {
        m.group(1): m.group(2) for m in _JS_VAR_RE.finditer(xml_text)
    }

    # Numeric fields
    def _js_float(key: str) -> float | None:
        raw = js_vars.get(key, "")
        # Strip trailing unit suffix (e.g. "165.0 cm" → "165.0")
        parts = raw.strip().split()
        return _to_float(parts[0]) if parts else None

    info.height_cm = _js_float("Height")
    info.weight_kg = _js_float("Weight")
    info.age = int(a) if (a := _to_float(js_vars.get("Age", ""))) is not None else None
    info.sex = js_vars.get("PatientSex") or None
    info.scan_date = js_vars.get("Scan") or None
    info.scan_protocol = js_vars.get("AnalProtocol") or None
    info.scan_mode = js_vars.get("ScanMode") or None

    # ACF / BCF from TableNote1
    note_raw = js_vars.get("TableNote1", "")
    note = _strip_html(note_raw)
    if note:
        info.calibration_note = note

    m_acf = _ACF_RE.search(note_raw)
    m_bcf = _BCF_RE.search(note_raw)
    if m_acf:
        info.acf = _to_float(m_acf.group(1))
    if m_bcf:
        info.bcf = _to_float(m_bcf.group(1))
