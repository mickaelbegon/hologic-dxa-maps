"""DICOM object classification for the Hologic DXA pipeline.

Classification is based solely on DICOM metadata (no pixel access required).
The result feeds the inventory and determines which processing path to take.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass

import pydicom
from pydicom.dataset import Dataset

from hologic_dxa.dicom.private_tags import PRDiagnosis, inspect_private_tags


class DicomClass(enum.Enum):
    """Coarse classification of a DICOM object."""

    SECONDARY_CAPTURE = "Secondary Capture"
    DXA_IMAGE = "DXA Image"
    STRUCTURED_REPORT = "Structured Report"
    PARAMETRIC_MAP = "Parametric Map"
    HOLOGIC_ARCHIVE = "Hologic Archive (P/R data)"
    UNKNOWN = "Unknown"


class QuantitativeEligibility(enum.Enum):
    """Whether this object can contribute calibrated quantitative data."""

    ELIGIBLE = "eligible"
    """Object contains or references calibrated quantitative values."""

    NOT_ELIGIBLE = "not_eligible"
    """Object is display-only or classification is unclear."""

    REQUIRES_CALIBRATION = "requires_calibration"
    """Raw signal present but calibration factors are missing."""

    UNKNOWN = "unknown"
    """Insufficient information to determine eligibility."""


# DICOM SOP Class UID constants
_SC_UID = "1.2.840.10008.5.1.4.1.1.7"
_SC_MULTI_UIDS = {
    "1.2.840.10008.5.1.4.1.1.7.1",
    "1.2.840.10008.5.1.4.1.1.7.2",
    "1.2.840.10008.5.1.4.1.1.7.3",
    "1.2.840.10008.5.1.4.1.1.7.4",
}
_SR_PREFIX = "1.2.840.10008.5.1.4.1.1.88."
_PARAMETRIC_MAP_UID = "1.2.840.10008.5.1.4.1.1.30"

# Hologic uses "OT" (Other) modality or "BMD" for DXA images
_DXA_MODALITIES = {"BMD", "OT", "DX"}
_DXA_MANUFACTURERS = {"hologic", "lunar", "norland"}


@dataclass(frozen=True)
class ClassificationResult:
    sop_class_uid: str
    modality: str
    manufacturer: str
    dicom_class: DicomClass
    pr_diagnosis: PRDiagnosis
    quantitative_eligibility: QuantitativeEligibility
    has_real_world_value_mapping: bool
    has_float_pixel_data: bool
    diagnostic: str


def classify(ds: Dataset) -> ClassificationResult:
    """Classify a pydicom Dataset without accessing pixel data.

    Parameters
    ----------
    ds:
        Dataset opened with stop_before_pixels=True (or full dataset).

    Returns
    -------
    ClassificationResult
    """
    sop_uid = str(getattr(ds, "SOPClassUID", ""))
    modality = str(getattr(ds, "Modality", "")).upper()
    manufacturer = str(getattr(ds, "Manufacturer", "")).lower()

    pr_result = inspect_private_tags(ds)

    # --- Classify DICOM class ---
    dicom_class = _classify_sop(sop_uid, modality, manufacturer, pr_result.diagnosis)

    # --- Check for quantitative metadata ---
    has_rwvm = (
        hasattr(ds, "RealWorldValueMappingSequence")
        and len(ds.RealWorldValueMappingSequence) > 0
    )
    has_float_pixels = (
        hasattr(ds, "FloatPixelData") or hasattr(ds, "DoubleFloatPixelData")
    )

    # --- Determine quantitative eligibility ---
    eligibility, diagnostic = _determine_eligibility(
        dicom_class, pr_result, has_rwvm, has_float_pixels, ds
    )

    return ClassificationResult(
        sop_class_uid=sop_uid,
        modality=modality,
        manufacturer=manufacturer,
        dicom_class=dicom_class,
        pr_diagnosis=pr_result.diagnosis,
        quantitative_eligibility=eligibility,
        has_real_world_value_mapping=has_rwvm,
        has_float_pixel_data=has_float_pixels,
        diagnostic=diagnostic,
    )


def _classify_sop(
    sop_uid: str,
    modality: str,
    manufacturer: str,
    pr_diagnosis: PRDiagnosis,
) -> DicomClass:
    if sop_uid == _PARAMETRIC_MAP_UID:
        return DicomClass.PARAMETRIC_MAP

    if sop_uid.startswith(_SR_PREFIX):
        return DicomClass.STRUCTURED_REPORT

    if pr_diagnosis == PRDiagnosis.PR_PRESENT:
        return DicomClass.HOLOGIC_ARCHIVE

    if sop_uid in (_SC_UID, *_SC_MULTI_UIDS):
        return DicomClass.SECONDARY_CAPTURE

    # Heuristic: DXA images often use OT or BMD modality
    if modality in _DXA_MODALITIES and any(m in manufacturer for m in _DXA_MANUFACTURERS):
        return DicomClass.DXA_IMAGE

    return DicomClass.UNKNOWN


def _determine_eligibility(
    dicom_class: DicomClass,
    pr_result: object,
    has_rwvm: bool,
    has_float_pixels: bool,
    ds: Dataset,
) -> tuple[QuantitativeEligibility, str]:
    if dicom_class == DicomClass.PARAMETRIC_MAP:
        if has_rwvm or has_float_pixels:
            return (
                QuantitativeEligibility.ELIGIBLE,
                "DICOM Parametric Map with quantitative pixel data.",
            )
        return (
            QuantitativeEligibility.NOT_ELIGIBLE,
            "Parametric Map missing RealWorldValueMappingSequence and float pixels.",
        )

    if dicom_class == DicomClass.STRUCTURED_REPORT:
        return (
            QuantitativeEligibility.ELIGIBLE,
            "Structured Report may contain regional quantitative results.",
        )

    if dicom_class == DicomClass.HOLOGIC_ARCHIVE:
        return (
            QuantitativeEligibility.REQUIRES_CALIBRATION,
            "Hologic archive P/R data present. Calibration required before "
            "quantitative reconstruction. See docs/limitations.md.",
        )

    if dicom_class == DicomClass.SECONDARY_CAPTURE:
        return (
            QuantitativeEligibility.NOT_ELIGIBLE,
            "Secondary Capture images MUST NOT be used as quantitative maps. "
            "Pixel values are display-scaled and carry no physical unit. "
            "See constraint #1 in docs/limitations.md.",
        )

    return (
        QuantitativeEligibility.UNKNOWN,
        f"Unrecognised DICOM class '{dicom_class.value}'. "
        "Manual inspection required.",
    )
