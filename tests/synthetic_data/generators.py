"""Synthetic DICOM dataset generators for testing.

All generated datasets use synthetic/random identifiers and contain no
protected health information (PHI). No real patient data is used anywhere
in this module.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

# ---------------------------------------------------------------------------
# SOP Class UID constants
# ---------------------------------------------------------------------------

_SC_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.7"
_SR_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.88.33"  # Comprehensive SR
_PM_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.30"      # Parametric Map


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_file_meta(sop_class_uid: str, sop_instance_uid: str) -> FileMetaDataset:
    meta = FileMetaDataset()
    meta.FileMetaInformationVersion = b"\x00\x01"
    meta.MediaStorageSOPClassUID = sop_class_uid
    meta.MediaStorageSOPInstanceUID = sop_instance_uid
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    return meta


def _base_dataset(sop_class_uid: str, modality: str = "OT") -> Dataset:
    """Minimal Dataset with standard DICOM identifiers — no PHI."""
    ds = Dataset()
    sop_uid = generate_uid()
    ds.SOPClassUID = sop_class_uid
    ds.SOPInstanceUID = sop_uid
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()
    ds.Modality = modality
    ds.Manufacturer = "SYNTHETIC"
    ds.ManufacturerModelName = "TestGenerator"
    ds.SoftwareVersions = "0.1.0"
    ds.file_meta = _make_file_meta(sop_class_uid, sop_uid)
    ds.preamble = b"\x00" * 128
    return ds


# ---------------------------------------------------------------------------
# Public generators
# ---------------------------------------------------------------------------

def make_secondary_capture(
    rows: int = 100,
    cols: int = 80,
    pixel_spacing_mm: tuple[float, float] = (1.0, 1.0),
) -> pydicom.Dataset:
    """Minimal Secondary Capture DICOM with synthetic pixel data."""
    ds = _base_dataset(_SC_SOP_CLASS, modality="OT")

    ds.Rows = rows
    ds.Columns = cols
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.PixelSpacing = list(pixel_spacing_mm)

    rng = np.random.default_rng(seed=42)
    pixel_array = rng.integers(0, 4096, size=(rows, cols), dtype=np.uint16)
    ds.PixelData = pixel_array.tobytes()

    return ds


def make_hologic_archive(
    rows: int = 100,
    cols: int = 80,
    p_data: bytes = b"FAKE_P_FILE_CONTENT",
    r_data: bytes = b"FAKE_R_FILE_CONTENT",
    p_declared_length: int | None = None,
    r_declared_length: int | None = None,
    add_padding_byte: bool = False,
    include_p: bool = True,
    include_r: bool = True,
) -> pydicom.Dataset:
    """DICOM with Hologic private P/R tags in group (0023).

    Parameters
    ----------
    p_declared_length:
        Explicit declared length stored in the P_FILE_LENGTH tag.
        ``None`` means use ``len(p_data)`` (before any padding).
    r_declared_length:
        Explicit declared length for the R file.
    add_padding_byte:
        When ``True``, one null byte is appended to each blob to simulate
        DICOM odd-length padding. The declared lengths are set to the
        *original* lengths so the padding byte triggers no warning.
    include_p:
        When ``False``, the P file tags are omitted entirely (simulates
        partial archive).
    include_r:
        When ``False``, the R file tags are omitted entirely (simulates
        partial archive).
    """
    ds = _base_dataset(_SC_SOP_CLASS, modality="BMD")
    ds.Manufacturer = "HOLOGIC, Inc."
    ds.ManufacturerModelName = "Discovery"

    ds.Rows = rows
    ds.Columns = cols

    # Compute actual blobs and declared lengths
    actual_p = (p_data + b"\x00") if add_padding_byte else p_data
    actual_r = (r_data + b"\x00") if add_padding_byte else r_data

    # Declared lengths refer to the original (non-padded) content length
    p_len = p_declared_length if p_declared_length is not None else len(p_data)
    r_len = r_declared_length if r_declared_length is not None else len(r_data)

    # Register private block — creator goes to (0023,0010)
    block = ds.private_block(0x0023, "HOLOGIC, Inc.", create=True)
    block.add_new(0x00, "LO", "HOLOGIC_ENCODING_V1")  # (0023,1000) ENCODING_SCHEME
    block.add_new(0x01, "LO", "scan.p")               # (0023,1001) P_FILE_NAME

    if include_p:
        block.add_new(0x02, "OB", actual_p)           # (0023,1002) P_FILE_DATA
        block.add_new(0x03, "UL", p_len)              # (0023,1003) P_FILE_LENGTH

    if include_r:
        block.add_new(0x04, "OB", actual_r)           # (0023,1004) R_FILE_DATA
        block.add_new(0x05, "UL", r_len)              # (0023,1005) R_FILE_LENGTH

    return ds


def make_structured_report(
    measurements: list[dict] | None = None,
) -> pydicom.Dataset:
    """Minimal DICOM SR with one NUM content item per measurement.

    Each measurement dict: {code_value, coding_scheme, code_meaning,
    value, unit_code, unit_meaning}.
    """
    ds = _base_dataset(_SR_SOP_CLASS, modality="SR")

    if measurements is None:
        measurements = [
            {
                "code_value": "122492",
                "coding_scheme": "DCM",
                "code_meaning": "Bone Mineral Density",
                "value": 1.23,
                "unit_code": "mg/cm2",
                "unit_meaning": "Milligrams per square centimeter",
            }
        ]

    content_items = []
    for meas in measurements:
        item = Dataset()
        item.RelationshipType = "CONTAINS"
        item.ValueType = "NUM"

        concept = Dataset()
        concept.CodeValue = meas["code_value"]
        concept.CodingSchemeDesignator = meas["coding_scheme"]
        concept.CodeMeaning = meas["code_meaning"]
        item.ConceptNameCodeSequence = Sequence([concept])

        mvsq = Dataset()
        mvsq.NumericValue = str(meas["value"])
        unit = Dataset()
        unit.CodeValue = meas["unit_code"]
        unit.CodingSchemeDesignator = "UCUM"
        unit.CodeMeaning = meas["unit_meaning"]
        mvsq.MeasurementUnitsCodeSequence = Sequence([unit])
        item.MeasuredValueSequence = Sequence([mvsq])

        content_items.append(item)

    ds.ContentSequence = Sequence(content_items)

    root_code = Dataset()
    root_code.CodeValue = "126000"
    root_code.CodingSchemeDesignator = "DCM"
    root_code.CodeMeaning = "Imaging Measurement Report"
    ds.ConceptNameCodeSequence = Sequence([root_code])
    ds.ValueType = "CONTAINER"
    ds.ContinuityOfContent = "SEPARATE"

    return ds


def make_parametric_map(
    rows: int = 50,
    cols: int = 40,
    float_values: np.ndarray | None = None,
    unit_code: str = "mg/cm2",
    pixel_spacing_mm: tuple[float, float] = (1.0, 1.0),
    include_rwvm: bool = True,
    include_float_pixels: bool = True,
) -> pydicom.Dataset:
    """DICOM Parametric Map with FloatPixelData and RealWorldValueMappingSequence.

    Parameters
    ----------
    include_rwvm:
        When ``False``, the RealWorldValueMappingSequence is omitted.
    include_float_pixels:
        When ``False``, FloatPixelData is omitted (creates a map with no
        quantitative pixel data — classified as NOT_ELIGIBLE).
    """
    ds = _base_dataset(_PM_SOP_CLASS, modality="OT")

    if float_values is None:
        rng = np.random.default_rng(seed=0)
        float_values = rng.uniform(0.0, 2.0, size=(rows, cols)).astype(np.float32)
    else:
        float_values = float_values.astype(np.float32)
        rows, cols = float_values.shape

    ds.Rows = rows
    ds.Columns = cols
    ds.PixelSpacing = list(pixel_spacing_mm)
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = 3
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"

    if include_float_pixels:
        ds.FloatPixelData = float_values.tobytes()

    if include_rwvm:
        rwvm_item = Dataset()
        rwvm_item.RealWorldValueFirstValueMapped = 0
        rwvm_item.RealWorldValueLastValueMapped = 65535
        rwvm_item.RealWorldValueSlope = 1.0
        rwvm_item.RealWorldValueIntercept = 0.0
        rwvm_item.LUTExplanation = f"Areal density ({unit_code})"
        rwvm_item.LUTLabel = "Areal density"

        unit_item = Dataset()
        unit_item.CodeValue = unit_code
        unit_item.CodingSchemeDesignator = "UCUM"
        unit_item.CodeMeaning = f"Areal density in {unit_code}"
        rwvm_item.MeasurementUnitsCodeSequence = Sequence([unit_item])

        ds.RealWorldValueMappingSequence = Sequence([rwvm_item])

    return ds


def make_truncated_pr_archive() -> pydicom.Dataset:
    """Archive DICOM where P file data is shorter than declared length.

    The declared length exceeds the actual length by 10 bytes, which is
    outside the ±1 tolerance for odd-length DICOM padding. This triggers
    a length-mismatch warning in ``inspect_private_tags``.
    """
    p_data = b"SHORT_P_DATA"
    return make_hologic_archive(
        p_data=p_data,
        r_data=b"VALID_R_DATA",
        p_declared_length=len(p_data) + 10,
    )


def make_empty_pr_archive() -> pydicom.Dataset:
    """Archive DICOM where P and R blobs are present but zero-length."""
    return make_hologic_archive(
        p_data=b"",
        r_data=b"",
    )


def save_to_tmpdir(
    ds: pydicom.Dataset,
    tmp_path: Path,
    filename: str,
) -> Path:
    """Save synthetic DICOM to a temp path. Returns the file path."""
    file_path = tmp_path / filename
    ds.save_as(str(file_path), write_like_original=False)
    return file_path
