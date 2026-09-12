"""Hologic APEX Structured Report parser.

Design principles
-----------------
* DICOM concept codes (CodeValue + CodingSchemeDesignator) are used as
  identifiers — never English label strings, which vary by locale and APEX version.
* The ContentSequence is traversed recursively.  Every NUM-type item is captured
  with its full code context, unit, anatomic region, and tree path.
* Output is long-format: one row per measurement, preserving all code fields so
  downstream code can filter by code rather than by string label.
* Non-SR datasets return an empty list rather than raising.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pydicom
import structlog

log = structlog.get_logger(__name__)

_SR_SOP_PREFIX = "1.2.840.10008.5.1.4.1.1.88."


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class SRMeasurement:
    """One numeric content item extracted from a DICOM SR."""

    code_value: str
    coding_scheme: str
    code_meaning: str
    value: float | None
    unit_code: str | None
    unit_meaning: str | None
    anatomic_region_code: str | None
    anatomic_region_meaning: str | None
    tracking_id: str | None
    path_in_tree: str


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_structured_report(ds: pydicom.Dataset) -> list[SRMeasurement]:
    """Recursively extract all NUM-type content items from a DICOM SR.

    Returns an empty list for non-SR datasets rather than raising, because
    callers typically process mixed DICOM directories.
    """
    sop_class = str(getattr(ds, "SOPClassUID", ""))
    if not sop_class.startswith(_SR_SOP_PREFIX):
        log.debug("not_a_structured_report", sop_class=sop_class)
        return []

    content_seq = getattr(ds, "ContentSequence", None)
    if content_seq is None:
        log.warning("sr_has_no_content_sequence")
        return []

    measurements: list[SRMeasurement] = []
    _traverse(content_seq, measurements, path="root", inherited_region=None)
    log.info("sr_parsed", measurement_count=len(measurements))
    return measurements


def sr_to_dataframe(measurements: list[SRMeasurement]) -> pd.DataFrame:
    """Convert a list of SRMeasurement to a long-format DataFrame.

    All code fields are preserved so downstream filters can match on
    code_value/coding_scheme rather than on human-readable strings.
    """
    if not measurements:
        return pd.DataFrame(
            columns=[
                "code_value", "coding_scheme", "code_meaning",
                "value", "unit_code", "unit_meaning",
                "anatomic_region_code", "anatomic_region_meaning",
                "tracking_id", "path_in_tree",
            ]
        )
    rows = [
        {
            "code_value": m.code_value,
            "coding_scheme": m.coding_scheme,
            "code_meaning": m.code_meaning,
            "value": m.value,
            "unit_code": m.unit_code,
            "unit_meaning": m.unit_meaning,
            "anatomic_region_code": m.anatomic_region_code,
            "anatomic_region_meaning": m.anatomic_region_meaning,
            "tracking_id": m.tracking_id,
            "path_in_tree": m.path_in_tree,
        }
        for m in measurements
    ]
    return pd.DataFrame(rows)


def load_sr_directory(sr_dir: Path) -> tuple[list[SRMeasurement], dict[str, Any]]:
    """Load all SR DICOMs from a directory.

    Returns
    -------
    measurements
        Flat list of all SRMeasurement objects from all files.
    tree
        Dict mapping filename → list of raw measurement dicts, for debugging.
    """
    all_measurements: list[SRMeasurement] = []
    tree: dict[str, Any] = {}

    dcm_files = sorted(sr_dir.glob("**/*.dcm"))
    if not dcm_files:
        dcm_files = sorted(sr_dir.glob("**/*"))

    for path in dcm_files:
        if not path.is_file():
            continue
        try:
            ds = pydicom.dcmread(str(path), stop_before_pixels=True)
        except Exception:  # noqa: BLE001
            log.debug("sr_dir_skip_unreadable")
            continue

        file_measurements = parse_structured_report(ds)
        if file_measurements:
            all_measurements.extend(file_measurements)
            tree[path.name] = [
                {"code_value": m.code_value, "code_meaning": m.code_meaning, "value": m.value}
                for m in file_measurements
            ]

    log.info(
        "sr_directory_loaded",
        file_count=len(tree),
        total_measurements=len(all_measurements),
    )
    return all_measurements, tree


# ---------------------------------------------------------------------------
# Recursive traversal
# ---------------------------------------------------------------------------

def _traverse(
    sequence: Any,
    out: list[SRMeasurement],
    path: str,
    inherited_region: tuple[str, str] | None,
) -> None:
    """Walk ContentSequence recursively, collecting NUM items."""
    for idx, item in enumerate(sequence):
        item_path = f"{path}/{idx}"
        value_type = str(getattr(item, "ValueType", "")).upper()

        region = _extract_anatomic_region(item) or inherited_region

        if value_type == "NUM":
            measurement = _extract_num_item(item, item_path, region)
            if measurement is not None:
                out.append(measurement)

        child_seq = getattr(item, "ContentSequence", None)
        if child_seq is not None:
            _traverse(child_seq, out, item_path, region)


def _extract_num_item(
    item: Any,
    path: str,
    region: tuple[str, str] | None,
) -> SRMeasurement | None:
    concept = _code_from_sequence(getattr(item, "ConceptNameCodeSequence", None))
    if concept is None:
        log.debug("num_item_missing_concept_code", path=path)
        return None

    numeric_value_seq = getattr(item, "MeasuredValueSequence", None)
    raw_value: float | None = None
    unit_code: str | None = None
    unit_meaning: str | None = None

    if numeric_value_seq and len(numeric_value_seq) > 0:
        mv = numeric_value_seq[0]
        try:
            raw_value = float(mv.NumericValue)
        except (AttributeError, ValueError, TypeError):
            pass
        unit = _code_from_sequence(getattr(mv, "MeasurementUnitsCodeSequence", None))
        if unit is not None:
            unit_code, _, unit_meaning = unit

    tracking_id = _safe_str_attr(item, "TrackingID")

    return SRMeasurement(
        code_value=concept[0],
        coding_scheme=concept[1],
        code_meaning=concept[2],
        value=raw_value,
        unit_code=unit_code,
        unit_meaning=unit_meaning,
        anatomic_region_code=region[0] if region else None,
        anatomic_region_meaning=region[1] if region else None,
        tracking_id=tracking_id,
        path_in_tree=path,
    )


def _extract_anatomic_region(item: Any) -> tuple[str, str] | None:
    region_seq = getattr(item, "AnatomicRegionSequence", None)
    if not region_seq or len(region_seq) == 0:
        return None
    code = _code_from_sequence(region_seq)
    if code is None:
        return None
    return (code[0], code[2])


def _code_from_sequence(seq: Any) -> tuple[str, str, str] | None:
    """Extract (CodeValue, CodingSchemeDesignator, CodeMeaning) from a code sequence."""
    if seq is None or len(seq) == 0:
        return None
    item = seq[0]
    code_value = _safe_str_attr(item, "CodeValue") or _safe_str_attr(item, "LongCodeValue")
    scheme = _safe_str_attr(item, "CodingSchemeDesignator")
    meaning = _safe_str_attr(item, "CodeMeaning")
    if not code_value or not scheme:
        return None
    return (code_value, scheme, meaning or "")


def _safe_str_attr(obj: Any, attr: str) -> str | None:
    val = getattr(obj, attr, None)
    if val is None:
        return None
    s = str(val).strip()
    return s if s else None
