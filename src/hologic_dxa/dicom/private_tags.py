"""Hologic-specific private DICOM tag definitions and inspection utilities.

The tag assignments below are based on empirical observation and publicly
available Hologic documentation fragments. They MUST be verified against
actual DICOM files or official Hologic documentation before relying on them
for data extraction.

Unknown group (0023,xxxx) tags are common in Hologic DXA DICOM archives.
Their presence does NOT imply that their content is decodable without
Hologic-provided documentation or SDK.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any

import pydicom
from pydicom.dataset import Dataset
from pydicom.tag import Tag


# ---------------------------------------------------------------------------
# Known Hologic private tag definitions
# ---------------------------------------------------------------------------

#: Private creator string that Hologic registers in group (0023,0010)
HOLOGIC_PRIVATE_CREATOR = "HOLOGIC, Inc."

#: Tags within private block (0023,10xx) — presence varies by APEX version
class HologicTag(enum.Enum):
    """Enumeration of known/suspected Hologic private tags in group 0023."""

    PRIVATE_CREATOR = Tag(0x0023, 0x0010)
    """Private creator identifier string."""

    ENCODING_SCHEME = Tag(0x0023, 0x1000)
    """Encoding scheme descriptor for the embedded P and R files."""

    P_FILE_NAME = Tag(0x0023, 0x1001)
    """Original name of the P (projection) file as stored by APEX."""

    P_FILE_DATA = Tag(0x0023, 0x1002)
    """Binary content of the P file."""

    P_FILE_LENGTH = Tag(0x0023, 0x1003)
    """Declared byte length of the P file content."""

    R_FILE_DATA = Tag(0x0023, 0x1004)
    """Binary content of the R file."""

    R_FILE_LENGTH = Tag(0x0023, 0x1005)
    """Declared byte length of the R file content."""


# ---------------------------------------------------------------------------
# Diagnosis result
# ---------------------------------------------------------------------------

class PRDiagnosis(enum.Enum):
    """Capability classification for Hologic P/R data in a DICOM object."""

    PR_PRESENT = "PR_PRESENT"
    """Both P and R file data tags found with non-zero content."""

    PR_PARTIAL = "PR_PARTIAL"
    """Some but not all expected P/R tags present."""

    PRIVATE_TAGS_PRESENT_BUT_UNKNOWN = "PRIVATE_TAGS_PRESENT_BUT_UNKNOWN"
    """Private group (0023) exists but tags do not match known Hologic schema."""

    ARCHIVE_DATA_ABSENT = "ARCHIVE_DATA_ABSENT"
    """No private tags in group (0023) detected."""

    SECONDARY_CAPTURE_ONLY = "SECONDARY_CAPTURE_ONLY"
    """Object is a Secondary Capture; no archive data present."""

    STRUCTURED_REPORT_ONLY = "STRUCTURED_REPORT_ONLY"
    """Object is a Structured Report; no image or archive data present."""

    PARAMETRIC_MAP_PRESENT = "PARAMETRIC_MAP_PRESENT"
    """Object is a DICOM Parametric Map with quantitative pixel data."""

    UNSUPPORTED = "UNSUPPORTED"
    """Object type is not recognised or cannot be processed."""


@dataclass
class PrivateTagInspectionResult:
    """Results of inspecting Hologic private tags in one DICOM object."""

    sop_instance_uid: str
    diagnosis: PRDiagnosis
    private_creator: str | None = None
    encoding_scheme: str | None = None
    p_file_name: str | None = None
    p_declared_length: int | None = None
    p_actual_length: int | None = None
    r_declared_length: int | None = None
    r_actual_length: int | None = None
    warnings: list[str] = field(default_factory=list)
    raw_private_tags: dict[str, Any] = field(default_factory=dict)


def inspect_private_tags(ds: Dataset) -> PrivateTagInspectionResult:
    """Inspect Hologic private tags in a pydicom Dataset.

    Does NOT decode the binary content of P/R files — that is the
    responsibility of :mod:`hologic_dxa.dicom.extract_pr`.

    Parameters
    ----------
    ds:
        A pydicom Dataset already loaded (pixels need not be present).

    Returns
    -------
    PrivateTagInspectionResult
        Inspection findings. Never raises; records anomalies as warnings.
    """
    sop_uid = str(getattr(ds, "SOPInstanceUID", "UNKNOWN"))
    warnings: list[str] = []
    raw: dict[str, Any] = {}

    # --- Private creator ---
    private_creator: str | None = None
    creator_tag = HologicTag.PRIVATE_CREATOR.value
    if creator_tag in ds:
        try:
            private_creator = str(ds[creator_tag].value)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"Could not read private creator: {exc}")

    # --- Check presence of group 0023 at all ---
    has_group_0023 = any(
        tag.group == 0x0023 for tag in ds.keys()
    )

    if not has_group_0023:
        return PrivateTagInspectionResult(
            sop_instance_uid=sop_uid,
            diagnosis=_diagnose_no_private(ds),
            warnings=warnings,
        )

    # --- Inventory all found Hologic tags ---
    found: set[HologicTag] = set()
    for ht in HologicTag:
        if ht.value in ds:
            found.add(ht)
            elem = ds[ht.value]
            # Store only metadata (length, VR) — not raw bytes — to avoid
            # accidental logging of large binary blobs
            raw[str(ht.name)] = {
                "VR": elem.VR,
                "length": len(elem.value) if hasattr(elem.value, "__len__") else None,
            }

    # Refuse to confirm Hologic creator if string does not match
    if private_creator is not None and HOLOGIC_PRIVATE_CREATOR not in private_creator:
        warnings.append(
            f"Private creator '{private_creator}' does not match expected "
            f"'{HOLOGIC_PRIVATE_CREATOR}'. Tag layout may differ."
        )

    encoding_scheme = _safe_str(ds, HologicTag.ENCODING_SCHEME.value)
    p_file_name = _safe_str(ds, HologicTag.P_FILE_NAME.value)

    p_declared = _safe_int(ds, HologicTag.P_FILE_LENGTH.value, warnings)
    p_actual = _blob_length(ds, HologicTag.P_FILE_DATA.value)

    r_declared = _safe_int(ds, HologicTag.R_FILE_LENGTH.value, warnings)
    r_actual = _blob_length(ds, HologicTag.R_FILE_DATA.value)

    # Length cross-checks
    if p_declared is not None and p_actual is not None and p_declared != p_actual:
        # DICOM odd-length padding: actual may be declared + 1
        if abs(p_declared - p_actual) != 1:
            warnings.append(
                f"P file length mismatch: declared={p_declared}, actual={p_actual}"
            )
    if r_declared is not None and r_actual is not None and r_declared != r_actual:
        if abs(r_declared - r_actual) != 1:
            warnings.append(
                f"R file length mismatch: declared={r_declared}, actual={r_actual}"
            )

    # Diagnose
    has_p = HologicTag.P_FILE_DATA in found and p_actual and p_actual > 0
    has_r = HologicTag.R_FILE_DATA in found and r_actual and r_actual > 0

    if has_p and has_r:
        diagnosis = PRDiagnosis.PR_PRESENT
    elif has_p or has_r:
        diagnosis = PRDiagnosis.PR_PARTIAL
        warnings.append("Only one of P or R file data tags found.")
    elif found:
        diagnosis = PRDiagnosis.PRIVATE_TAGS_PRESENT_BUT_UNKNOWN
        warnings.append(
            "Private group 0023 present but does not contain expected P/R data tags. "
            "Possibly a different APEX version or a non-archive object."
        )
    else:
        diagnosis = PRDiagnosis.ARCHIVE_DATA_ABSENT

    return PrivateTagInspectionResult(
        sop_instance_uid=sop_uid,
        diagnosis=diagnosis,
        private_creator=private_creator,
        encoding_scheme=encoding_scheme,
        p_file_name=p_file_name,
        p_declared_length=p_declared,
        p_actual_length=p_actual,
        r_declared_length=r_declared,
        r_actual_length=r_actual,
        warnings=warnings,
        raw_private_tags=raw,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _diagnose_no_private(ds: Dataset) -> PRDiagnosis:
    sop_class = str(getattr(ds, "SOPClassUID", ""))
    # DICOM SOP Class UIDs (subset)
    _SECONDARY_CAPTURE = "1.2.840.10008.5.1.4.1.1.7"
    _SC_MULTIFRAME = "1.2.840.10008.5.1.4.1.1.7.2"
    _SR_ENHANCED = "1.2.840.10008.5.1.4.1.1.88"
    _SR_PREFIX = "1.2.840.10008.5.1.4.1.1.88."
    _PARAMETRIC_MAP = "1.2.840.10008.5.1.4.1.1.30"

    if sop_class in (_SECONDARY_CAPTURE, _SC_MULTIFRAME):
        return PRDiagnosis.SECONDARY_CAPTURE_ONLY
    if sop_class.startswith(_SR_PREFIX) or sop_class == _SR_ENHANCED:
        return PRDiagnosis.STRUCTURED_REPORT_ONLY
    if sop_class == _PARAMETRIC_MAP:
        return PRDiagnosis.PARAMETRIC_MAP_PRESENT
    return PRDiagnosis.ARCHIVE_DATA_ABSENT


def _safe_str(ds: Dataset, tag: Tag) -> str | None:
    if tag not in ds:
        return None
    try:
        return str(ds[tag].value)
    except Exception:  # noqa: BLE001
        return None


def _safe_int(ds: Dataset, tag: Tag, warnings: list[str]) -> int | None:
    if tag not in ds:
        return None
    try:
        val = ds[tag].value
        return int(val)
    except (ValueError, TypeError) as exc:
        warnings.append(f"Could not parse tag {tag} as int: {exc}")
        return None


def _blob_length(ds: Dataset, tag: Tag) -> int | None:
    if tag not in ds:
        return None
    try:
        val = ds[tag].value
        return len(val) if isinstance(val, (bytes, bytearray)) else None
    except Exception:  # noqa: BLE001
        return None
