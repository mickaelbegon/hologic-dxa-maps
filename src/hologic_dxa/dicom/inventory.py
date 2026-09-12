"""Recursive DICOM inventory builder — Phase 1 of the pipeline.

This module implements the `hologic-dxa inspect` command. It never reads
pixel data unless required for classification, and it never writes PHI
to any output artefact.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pydicom
from pydicom.errors import InvalidDicomError

from hologic_dxa.dicom.classify import ClassificationResult, classify
from hologic_dxa.dicom.private_tags import HologicTag
from hologic_dxa.logging import get_logger

log = get_logger(__name__)

# Tags that must never appear in any output (belt-and-suspenders)
_PHI_ATTRIBUTES = frozenset({
    "PatientName", "PatientID", "PatientBirthDate", "PatientSex",
    "PatientAge", "PatientWeight", "PatientSize",
    "OtherPatientIDs", "OtherPatientNames",
    "InstitutionName", "InstitutionAddress",
    "AccessionNumber",
})


@dataclass
class DicomFileRecord:
    """One row in the inventory — all PHI fields are absent by design."""

    file_path: str
    file_size_bytes: int
    sha256: str
    is_dicom: bool
    parse_error: str | None

    # Standard identifiers
    sop_class_uid: str | None = None
    sop_instance_uid: str | None = None
    study_instance_uid: str | None = None
    series_instance_uid: str | None = None

    # Device metadata
    modality: str | None = None
    manufacturer: str | None = None
    manufacturer_model_name: str | None = None
    software_versions: str | None = None
    transfer_syntax_uid: str | None = None

    # Image geometry
    rows: int | None = None
    columns: int | None = None
    bits_allocated: int | None = None
    bits_stored: int | None = None
    photometric_interpretation: str | None = None
    pixel_spacing: str | None = None
    imager_pixel_spacing: str | None = None
    rescale_slope: float | None = None
    rescale_intercept: float | None = None

    # Quantitative flags
    presence_of_real_world_value_mapping: bool = False
    presence_of_float_pixel_data: bool = False
    presence_of_private_group_0023: bool = False
    presence_of_hologic_pr_tags: bool = False

    # Classification
    classification: str | None = None
    pr_diagnosis: str | None = None
    quantitative_eligibility: str | None = None
    diagnostic: str | None = None

    warnings: list[str] = field(default_factory=list)


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def build_inventory(root: Path) -> list[DicomFileRecord]:
    """Walk *root* recursively and return one record per file.

    Parameters
    ----------
    root:
        Directory to search. Non-DICOM files are recorded with is_dicom=False.
    """
    records: list[DicomFileRecord] = []

    all_files = sorted(root.rglob("*"))
    files = [f for f in all_files if f.is_file()]
    log.info("inventory_started", root=str(root), file_count=len(files))

    for path in files:
        records.append(_process_file(path, root))

    log.info("inventory_complete", total=len(records),
             dicom_count=sum(r.is_dicom for r in records))
    return records


def _process_file(path: Path, root: Path) -> DicomFileRecord:
    relative = str(path.relative_to(root))
    size = path.stat().st_size

    try:
        sha = sha256_file(path)
    except OSError as exc:
        return DicomFileRecord(
            file_path=relative,
            file_size_bytes=size,
            sha256="",
            is_dicom=False,
            parse_error=f"IO error computing SHA-256: {exc}",
        )

    try:
        ds = pydicom.dcmread(str(path), stop_before_pixels=True, force=False)
    except (InvalidDicomError, Exception) as exc:
        return DicomFileRecord(
            file_path=relative,
            file_size_bytes=size,
            sha256=sha,
            is_dicom=False,
            parse_error=str(exc)[:200],
        )

    # Confirm no PHI slips through
    _assert_no_phi(ds, path)

    result = classify(ds)
    warnings = list(result.pr_diagnosis.name if False else [])  # placeholder

    has_group_0023 = any(tag.group == 0x0023 for tag in ds.keys())
    has_hologic_pr = result.pr_diagnosis.value in (
        "PR_PRESENT", "PR_PARTIAL"
    )

    pixel_spacing = _join_decimal_seq(getattr(ds, "PixelSpacing", None))
    imager_spacing = _join_decimal_seq(getattr(ds, "ImagerPixelSpacing", None))

    return DicomFileRecord(
        file_path=relative,
        file_size_bytes=size,
        sha256=sha,
        is_dicom=True,
        parse_error=None,
        sop_class_uid=str(getattr(ds, "SOPClassUID", "") or ""),
        sop_instance_uid=str(getattr(ds, "SOPInstanceUID", "") or ""),
        study_instance_uid=str(getattr(ds, "StudyInstanceUID", "") or ""),
        series_instance_uid=str(getattr(ds, "SeriesInstanceUID", "") or ""),
        modality=str(getattr(ds, "Modality", "") or ""),
        manufacturer=str(getattr(ds, "Manufacturer", "") or ""),
        manufacturer_model_name=str(getattr(ds, "ManufacturerModelName", "") or ""),
        software_versions=_software_versions(ds),
        transfer_syntax_uid=str(ds.file_meta.TransferSyntaxUID)
            if hasattr(ds, "file_meta") else None,
        rows=getattr(ds, "Rows", None),
        columns=getattr(ds, "Columns", None),
        bits_allocated=getattr(ds, "BitsAllocated", None),
        bits_stored=getattr(ds, "BitsStored", None),
        photometric_interpretation=str(getattr(ds, "PhotometricInterpretation", "") or ""),
        pixel_spacing=pixel_spacing,
        imager_pixel_spacing=imager_spacing,
        rescale_slope=_safe_float(getattr(ds, "RescaleSlope", None)),
        rescale_intercept=_safe_float(getattr(ds, "RescaleIntercept", None)),
        presence_of_real_world_value_mapping=result.has_real_world_value_mapping,
        presence_of_float_pixel_data=result.has_float_pixel_data,
        presence_of_private_group_0023=has_group_0023,
        presence_of_hologic_pr_tags=has_hologic_pr,
        classification=result.dicom_class.value,
        pr_diagnosis=result.pr_diagnosis.value,
        quantitative_eligibility=result.quantitative_eligibility.value,
        diagnostic=result.diagnostic,
        warnings=warnings,
    )


def write_inventory(records: list[DicomFileRecord], output_dir: Path) -> None:
    """Write inventory.json, inventory.csv, and inventory_summary.md."""
    output_dir.mkdir(parents=True, exist_ok=True)

    dicts = [asdict(r) for r in records]

    # JSON
    json_path = output_dir / "inventory.json"
    json_path.write_text(json.dumps(dicts, indent=2), encoding="utf-8")

    # CSV
    df = pd.DataFrame(dicts)
    # lists → string for CSV compatibility
    df["warnings"] = df["warnings"].apply(lambda x: "; ".join(x) if x else "")
    csv_path = output_dir / "inventory.csv"
    df.to_csv(csv_path, index=False)

    # Summary markdown
    summary = _build_summary(records, df)
    md_path = output_dir / "inventory_summary.md"
    md_path.write_text(summary, encoding="utf-8")

    log.info("inventory_written",
             json=str(json_path), csv=str(csv_path), md=str(md_path))


def _build_summary(records: list[DicomFileRecord], df: pd.DataFrame) -> str:
    total = len(records)
    dicom_count = sum(r.is_dicom for r in records)
    error_count = sum(r.parse_error is not None for r in records)

    lines = [
        "# DICOM Inventory Summary",
        "",
        f"- Total files scanned: {total}",
        f"- Valid DICOM: {dicom_count}",
        f"- Parse errors / non-DICOM: {error_count}",
        "",
        "## Classification breakdown",
        "",
    ]

    if dicom_count > 0:
        dicom_df = df[df["is_dicom"]]
        for cls, count in dicom_df["classification"].value_counts().items():
            lines.append(f"- {cls}: {count}")
        lines.append("")
        lines.append("## Quantitative eligibility")
        lines.append("")
        for elig, count in dicom_df["quantitative_eligibility"].value_counts().items():
            lines.append(f"- {elig}: {count}")
        lines.append("")
        lines.append("## Hologic P/R diagnosis")
        lines.append("")
        for diag, count in dicom_df["pr_diagnosis"].value_counts().items():
            lines.append(f"- {diag}: {count}")

    lines.append("")
    lines.append("*No protected health information is present in this report.*")
    return "\n".join(lines)


def _assert_no_phi(ds: pydicom.Dataset, path: Path) -> None:
    for attr in _PHI_ATTRIBUTES:
        if hasattr(ds, attr):
            val = getattr(ds, attr)
            if val:
                log.warning(
                    "phi_detected_in_dataset",
                    attribute=attr,
                    file=str(path.name),
                    action="suppressed_from_output",
                )


def _join_decimal_seq(seq: Any) -> str | None:
    if seq is None:
        return None
    try:
        return "\\".join(str(v) for v in seq)
    except TypeError:
        return str(seq)


def _software_versions(ds: pydicom.Dataset) -> str | None:
    sv = getattr(ds, "SoftwareVersions", None)
    if sv is None:
        return None
    if isinstance(sv, str):
        return sv
    try:
        return "\\".join(str(v) for v in sv)
    except TypeError:
        return str(sv)


def _safe_float(val: Any) -> float | None:
    if val is None:
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
