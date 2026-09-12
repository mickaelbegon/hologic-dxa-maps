"""Export quantitative maps as DICOM Parametric Map objects using highdicom.

This export is ONLY activated when:
1. Units are known and expressed in g/cm²
2. Pixel geometry is defined (pixel spacing in mm)
3. Source DICOM UIDs are available for reference
4. Values have been validated against APEX regional results

Do not use this exporter for unvalidated experimental reconstructions.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from hologic_dxa.logging import get_logger
from hologic_dxa.providers.base import QuantitativeDataUnavailableError

if TYPE_CHECKING:
    from hologic_dxa.models import QuantitativeMap

log = get_logger(__name__)


def export_parametric_map_dicom(
    qmap: "QuantitativeMap",
    output_path: Path,
    source_sop_instance_uid: str,
    pixel_spacing_mm: tuple[float, float],
    validated: bool = False,
) -> None:
    """Export a QuantitativeMap as a DICOM Parametric Map.

    Parameters
    ----------
    validated:
        Must be explicitly set to True after mass conservation validation.
        Raises ValueError if False.

    Raises
    ------
    QuantitativeDataUnavailableError
        If units, geometry, or source UIDs are missing.
    ValueError
        If validated=False (prevents exporting unvalidated maps as DICOM).
    ImportError
        If highdicom is not installed.
    """
    if not validated:
        raise ValueError(
            "Cannot export unvalidated maps as DICOM Parametric Map. "
            "Run validate-maps first and set validated=True explicitly "
            "after reviewing the validation report."
        )

    if qmap.unit not in ("g/cm2", "g/cm²"):
        raise QuantitativeDataUnavailableError(
            f"DICOM PM export requires unit 'g/cm2', got '{qmap.unit}'. "
            "Convert to g/cm2 before exporting."
        )

    if pixel_spacing_mm[0] <= 0 or pixel_spacing_mm[1] <= 0:
        raise QuantitativeDataUnavailableError(
            "Pixel spacing must be positive. "
            "Cannot create a DICOM Parametric Map without valid geometry."
        )

    try:
        import highdicom as hd  # noqa: F401
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "DICOM PM export requires highdicom. "
            "Install with: pip install highdicom"
        ) from exc

    # Full highdicom PM construction is deferred until the highdicom API
    # is confirmed compatible with the project's DICOM requirements.
    # See: https://highdicom.readthedocs.io/en/latest/
    raise NotImplementedError(
        "DICOM Parametric Map export via highdicom is not yet implemented. "
        "The interface is defined and validation gates are in place. "
        "Implementation requires confirming: "
        "(1) highdicom ParametricMap API stability, "
        "(2) correct RWVM unit encoding for g/cm2, "
        "(3) source DICOM series reference format. "
        "Interim export: use export-hdf5 instead."
    )
