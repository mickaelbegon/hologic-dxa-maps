"""DICOM de-identification following DICOM PS 3.15 Annex E.

IMPORTANT: Always extract P/R files (dicom.extract_pr) BEFORE de-identifying.
De-identification may remove the private tags containing the P/R data.

De-identified copies are always written to a separate output directory.
Source files are never modified.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pydicom
from pydicom.dataset import Dataset
from pydicom.uid import generate_uid

from hologic_dxa.logging import get_logger
from hologic_dxa.provenance import PseudonymRegistry

log = get_logger(__name__)

# DICOM PS 3.15 Annex E — Basic Application Level Confidentiality Profile
# Attributes to remove or empty (subset — not exhaustive)
_REMOVE_ATTRIBUTES = [
    "PatientName",
    "PatientID",
    "PatientBirthDate",
    "PatientSex",
    "PatientAge",
    "PatientWeight",
    "PatientSize",
    "OtherPatientIDs",
    "OtherPatientNames",
    "OtherPatientIDsSequence",
    "InstitutionName",
    "InstitutionAddress",
    "InstitutionalDepartmentName",
    "ReferringPhysicianName",
    "PerformingPhysicianName",
    "RequestingPhysician",
    "OperatorsName",
    "AccessionNumber",
    "StudyID",
    "RequestedProcedureID",
    "ScheduledProcedureStepID",
]

# UIDs to replace with generated equivalents (maintain consistency within study)
_REPLACE_UIDS = [
    "StudyInstanceUID",
    "SeriesInstanceUID",
    "SOPInstanceUID",
    "FrameOfReferenceUID",
]


def deidentify_dataset(
    ds: Dataset,
    registry: PseudonymRegistry,
    uid_map: dict[str, str] | None = None,
) -> tuple[Dataset, dict[str, str]]:
    """Return a de-identified copy of the dataset.

    Parameters
    ----------
    ds:
        Original pydicom Dataset (not modified).
    registry:
        PseudonymRegistry for patient → research ID mapping.
    uid_map:
        Optional dict mapping original UID → replacement UID.
        Pass the same dict across calls to maintain UID consistency within a study.
        Modified in-place; returned for chaining.

    Returns
    -------
    (deidentified_ds, uid_map)
    """
    if uid_map is None:
        uid_map = {}

    ds_copy = ds.copy()

    # Replace patient ID with research pseudonym
    original_patient_id = str(getattr(ds_copy, "PatientID", "UNKNOWN"))
    research_id = registry.get_or_create(original_patient_id)
    log.info("deidentify_patient", research_id=research_id)

    # Remove PHI attributes
    for attr in _REMOVE_ATTRIBUTES:
        if hasattr(ds_copy, attr):
            delattr(ds_copy, attr)

    # Set research pseudonym as new PatientID (research use only)
    ds_copy.PatientID = research_id
    ds_copy.PatientName = research_id

    # Replace UIDs consistently
    for uid_attr in _REPLACE_UIDS:
        original = str(getattr(ds_copy, uid_attr, ""))
        if not original:
            continue
        if original not in uid_map:
            uid_map[original] = generate_uid()
        setattr(ds_copy, uid_attr, uid_map[original])

    # Remove private tags (group 0023 Hologic data) — caller must ensure
    # P/R extraction was done first
    ds_copy.remove_private_tags()
    log.info(
        "deidentify_complete",
        research_id=research_id,
        note="Private tags removed. Ensure PR extraction was done before this step.",
    )

    return ds_copy, uid_map


def deidentify_file(
    source_path: Path,
    output_dir: Path,
    registry: PseudonymRegistry,
    uid_map: dict[str, str] | None = None,
) -> tuple[Path, dict[str, str]]:
    """De-identify a single DICOM file and write to output_dir."""
    ds = pydicom.dcmread(str(source_path))
    ds_deident, uid_map_out = deidentify_dataset(ds, registry, uid_map)

    # Use a hash of the original SOP UID for the output filename — never
    # use the original patient-identifiable filename
    sop_uid = str(getattr(ds, "SOPInstanceUID", "UNKNOWN"))
    safe_name = hashlib.sha256(sop_uid.encode()).hexdigest()[:16] + ".dcm"
    output_path = output_dir / safe_name
    output_dir.mkdir(parents=True, exist_ok=True)
    pydicom.dcmwrite(str(output_path), ds_deident)
    log.info("deidentify_file_written", output=str(output_path))
    return output_path, uid_map_out
