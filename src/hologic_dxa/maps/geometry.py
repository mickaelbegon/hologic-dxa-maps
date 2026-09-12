"""Pixel coordinate and area calculations for DXA fan-beam scanners.

Critical limitation: In fan-beam DXA, pixel area is NOT constant — it varies
with distance from the beam centre. Without scanner geometry parameters (focal-spot
position, source-to-detector distance), only an approximate scalar pixel area
derived from PixelSpacing or ImagerPixelSpacing is available.  All functions in
this module operate under the flat-detector approximation unless otherwise stated.
"""

from __future__ import annotations

import numpy as np

from hologic_dxa.logging import get_logger

log = get_logger(__name__)


def pixel_area_cm2_from_spacing(
    row_spacing_mm: float, col_spacing_mm: float
) -> float:
    """Convert pixel spacing (mm) to pixel area (cm²). Flat-detector approximation.

    Parameters
    ----------
    row_spacing_mm:
        Centre-to-centre distance between adjacent rows in mm
        (first element of DICOM PixelSpacing or ImagerPixelSpacing).
    col_spacing_mm:
        Centre-to-centre distance between adjacent columns in mm.

    Returns
    -------
    float
        Pixel area in cm².

    Raises
    ------
    ValueError
        If either spacing value is not strictly positive.
    """
    if row_spacing_mm <= 0 or col_spacing_mm <= 0:
        raise ValueError(
            f"Pixel spacing values must be strictly positive; "
            f"got row_spacing_mm={row_spacing_mm}, col_spacing_mm={col_spacing_mm}. "
            "Verify that PixelSpacing or ImagerPixelSpacing was read correctly from "
            "the DICOM header."
        )
    area: float = (row_spacing_mm / 10.0) * (col_spacing_mm / 10.0)
    log.debug(
        "pixel_area_from_spacing",
        row_spacing_mm=row_spacing_mm,
        col_spacing_mm=col_spacing_mm,
        pixel_area_cm2=area,
    )
    return area


def pixel_coordinates_mm(
    rows: int,
    cols: int,
    row_spacing_mm: float,
    col_spacing_mm: float,
    origin_mm: tuple[float, float] = (0.0, 0.0),
) -> tuple[np.ndarray, np.ndarray]:
    """Return (x_coords, y_coords) arrays of shape (rows, cols) in mm.

    x varies along columns (lateral), y along rows (superior-inferior).

    Parameters
    ----------
    rows, cols:
        Image dimensions in pixels.
    row_spacing_mm, col_spacing_mm:
        Pixel spacing in mm.
    origin_mm:
        (x, y) position of the top-left pixel centre in mm.
        Defaults to (0, 0).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        x_coords and y_coords, each of shape (rows, cols) in mm.
        x_coords[r, c] is the lateral coordinate of pixel (r, c);
        y_coords[r, c] is the superior-inferior coordinate.
    """
    x_1d = origin_mm[0] + np.arange(cols, dtype=np.float64) * col_spacing_mm
    y_1d = origin_mm[1] + np.arange(rows, dtype=np.float64) * row_spacing_mm
    x_coords, y_coords = np.meshgrid(x_1d, y_1d)
    return x_coords, y_coords


def validate_pixel_area(pixel_area: np.ndarray | float | None) -> None:
    """Raise ValueError if pixel_area is None, zero, or negative.

    Parameters
    ----------
    pixel_area:
        Scalar or array of pixel areas in cm².

    Raises
    ------
    ValueError
        If pixel_area is None, empty, contains zeros, or contains negative values.
        The message explains what is needed to compute a valid pixel area.
    """
    if pixel_area is None:
        raise ValueError(
            "pixel_area is None. Supply PixelSpacing or ImagerPixelSpacing from the "
            "DICOM header and call pixel_area_cm2_from_spacing() to obtain a "
            "flat-detector approximation, or provide scanner geometry parameters "
            "(focal-spot position, source-to-detector distance) for an exact "
            "fan-beam pixel area."
        )
    arr = np.asarray(pixel_area, dtype=np.float64)
    if arr.size == 0:
        raise ValueError("pixel_area is an empty array.")
    min_val = float(np.nanmin(arr))
    if min_val <= 0.0:
        raise ValueError(
            f"pixel_area must be strictly positive in all elements; "
            f"got minimum value {min_val:.6g} cm². "
            "Check that PixelSpacing or ImagerPixelSpacing was parsed correctly "
            "and contains values in millimetres."
        )
