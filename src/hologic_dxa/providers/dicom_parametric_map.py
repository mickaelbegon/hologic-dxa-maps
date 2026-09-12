"""DICOM Parametric Map provider.

Reads quantitative pixel maps from DICOM Parametric Map objects. This is
the highest-confidence source of calibrated data because the DICOM standard
requires explicit unit encoding via RealWorldValueMappingSequence.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom

from hologic_dxa.dicom.parametric_map import load_parametric_map
from hologic_dxa.logging import get_logger
from hologic_dxa.maps.units import convert_to_g_cm2, validate_areal_density_unit
from hologic_dxa.providers.base import QuantitativeDataUnavailableError

log = get_logger(__name__)

_PARAMETRIC_MAP_UID = "1.2.840.10008.5.1.4.1.1.30"


class DicomParametricMapProvider:
    """Provider for DICOM Parametric Map objects (SOPClassUID .30)."""

    name = "dicom_parametric_map"
    experimental = False

    def can_read(self, source: Path) -> bool:
        if not source.is_file():
            return False
        try:
            ds = pydicom.dcmread(str(source), stop_before_pixels=True)
            return str(getattr(ds, "SOPClassUID", "")) == _PARAMETRIC_MAP_UID
        except Exception:  # noqa: BLE001
            return False

    def load(self, source: Path) -> object:
        """Load a single DICOM Parametric Map.

        Returns a dict (not QuantitativeMapBundle) because a single PM
        contains only one tissue component. Callers must assemble a bundle
        from multiple PM files (fat, lean, bmc).

        Raises
        ------
        QuantitativeDataUnavailableError
            If unit, geometry, or pixel data is missing or unrecognised.
        """
        values, unit, metadata = load_parametric_map(source)

        try:
            validate_areal_density_unit(unit)
        except ValueError as exc:
            raise QuantitativeDataUnavailableError(
                f"DICOM PM unit '{unit}' is not a recognised areal density unit. "
                "Cannot load as a quantitative map. Details: " + str(exc)
            ) from exc

        values_g_cm2 = convert_to_g_cm2(values, unit)

        pixel_spacing = metadata.get("pixel_spacing_mm")
        if not pixel_spacing or pixel_spacing[0] is None:
            raise QuantitativeDataUnavailableError(
                "Parametric Map has no PixelSpacing tag. "
                "Cannot compute pixel area for mass integration."
            )

        log.info(
            "dicom_pm_loaded",
            rows=metadata["rows"],
            cols=metadata["columns"],
            unit=unit,
            source=source.name,
        )

        return {
            "values_g_cm2": values_g_cm2,
            "unit": "g/cm2",
            "metadata": metadata,
            "pixel_spacing_mm": pixel_spacing,
        }

    def describe_capabilities(self) -> dict:
        return {
            "status": "implemented",
            "experimental": False,
            "supports": [
                "FloatPixelData",
                "DoubleFloatPixelData",
                "RealWorldValueMappingSequence",
            ],
            "limitations": [
                "Each PM file contains only one tissue component.",
                "Requires separate files for fat, lean, and BMC.",
                "Pixel area assumed constant (flat-detector approximation).",
            ],
        }
