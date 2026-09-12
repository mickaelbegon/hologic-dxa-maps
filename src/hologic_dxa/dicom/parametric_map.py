"""DICOM Parametric Map reader.

Reads quantitative pixel data from DICOM Parametric Map objects
(SOPClassUID 1.2.840.10008.5.1.4.1.1.30).

A Parametric Map is considered quantitative ONLY when:
- FloatPixelData or DoubleFloatPixelData is present
- RealWorldValueMappingSequence defines the unit relationship

Objects with neither of these are classified as NOT_ELIGIBLE and will
not be converted to QuantitativeMap objects.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import Dataset

from hologic_dxa.logging import get_logger
from hologic_dxa.providers.base import QuantitativeDataUnavailableError

log = get_logger(__name__)

_PARAMETRIC_MAP_UID = "1.2.840.10008.5.1.4.1.1.30"


def read_float_pixels(ds: Dataset) -> np.ndarray:
    """Extract float or double pixel data from a Parametric Map.

    Raises
    ------
    QuantitativeDataUnavailableError
        If neither FloatPixelData nor DoubleFloatPixelData is present.
    """
    if hasattr(ds, "FloatPixelData"):
        raw = ds.FloatPixelData
        dtype = np.float32
    elif hasattr(ds, "DoubleFloatPixelData"):
        raw = ds.DoubleFloatPixelData
        dtype = np.float64
    else:
        raise QuantitativeDataUnavailableError(
            "Parametric Map has neither FloatPixelData nor DoubleFloatPixelData. "
            "This object does not contain quantitative pixel data."
        )

    arr = np.frombuffer(raw, dtype=dtype).reshape(ds.Rows, ds.Columns)
    return arr.astype(np.float64)


def read_real_world_value_mapping(ds: Dataset) -> tuple[float, float, str]:
    """Extract slope, intercept, and unit from RealWorldValueMappingSequence.

    Returns
    -------
    (slope, intercept, unit_meaning)
        stored_value × slope + intercept = physical_value [unit_meaning]

    Raises
    ------
    QuantitativeDataUnavailableError
        If RWVM sequence is absent or malformed.
    """
    if not hasattr(ds, "RealWorldValueMappingSequence"):
        raise QuantitativeDataUnavailableError(
            "No RealWorldValueMappingSequence found. "
            "Cannot determine the physical unit of the pixel values. "
            "This object is not suitable as a quantitative map."
        )

    rwvm = ds.RealWorldValueMappingSequence[0]
    slope = float(getattr(rwvm, "RealWorldValueSlope", 1.0))
    intercept = float(getattr(rwvm, "RealWorldValueIntercept", 0.0))

    unit_meaning = "unknown"
    if hasattr(rwvm, "MeasurementUnitsCodeSequence"):
        unit_seq = rwvm.MeasurementUnitsCodeSequence[0]
        unit_meaning = str(getattr(unit_seq, "CodeMeaning", "unknown"))

    log.info(
        "rwvm_read",
        slope=slope,
        intercept=intercept,
        unit=unit_meaning,
    )
    return slope, intercept, unit_meaning


def load_parametric_map(path: Path) -> tuple[np.ndarray, str, dict]:
    """Load a DICOM Parametric Map and return (values_g_cm2, unit, metadata).

    Returns
    -------
    values:
        2D float64 array in physical units (after slope/intercept application).
    unit:
        Unit string from RWVM (e.g. "mg/cm2"). Caller must validate.
    metadata:
        Dict with SOPInstanceUID, PixelSpacing, etc.
    """
    ds = pydicom.dcmread(str(path))

    if str(getattr(ds, "SOPClassUID", "")) != _PARAMETRIC_MAP_UID:
        raise QuantitativeDataUnavailableError(
            f"'{path.name}' is not a DICOM Parametric Map "
            f"(SOPClassUID = {getattr(ds, 'SOPClassUID', 'missing')})."
        )

    pixels = read_float_pixels(ds)
    slope, intercept, unit = read_real_world_value_mapping(ds)
    physical_values = pixels * slope + intercept

    pixel_spacing = list(getattr(ds, "PixelSpacing", [None, None]))

    metadata = {
        "sop_instance_uid": str(getattr(ds, "SOPInstanceUID", "")),
        "study_instance_uid": str(getattr(ds, "StudyInstanceUID", "")),
        "series_instance_uid": str(getattr(ds, "SeriesInstanceUID", "")),
        "rows": ds.Rows,
        "columns": ds.Columns,
        "pixel_spacing_mm": pixel_spacing,
        "rwvm_slope": slope,
        "rwvm_intercept": intercept,
        "unit_from_rwvm": unit,
        "manufacturer": str(getattr(ds, "Manufacturer", "")),
        "software_versions": str(getattr(ds, "SoftwareVersions", "")),
    }

    return physical_values, unit, metadata
