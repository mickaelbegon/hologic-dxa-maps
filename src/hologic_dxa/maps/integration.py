"""Mass and inertia integration over quantitative DXA pixel maps.

All integrals are discrete sums:

    mass = Σ_i  ρ_i · A_i

where ρ_i is the areal density (g/cm²) and A_i is the pixel area (cm²).
Results are in grams; moments of inertia are in g·mm² (biomechanical convention).

NaN pixels within the mask are excluded via ``np.nansum`` so that partially
masked or edge-corrected images do not silently corrupt aggregate quantities.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from hologic_dxa.logging import get_logger
from hologic_dxa.maps.geometry import validate_pixel_area

if TYPE_CHECKING:
    from hologic_dxa.models import QuantitativeMapBundle

log = get_logger(__name__)


def regional_mass_g(
    areal_density: np.ndarray,
    pixel_area_cm2: np.ndarray | float,
    mask: np.ndarray,
) -> float:
    """Sum of areal_density * pixel_area over masked pixels.  Returns grams.

    Parameters
    ----------
    areal_density:
        2-D array of areal density values in g/cm².
    pixel_area_cm2:
        Scalar or array (broadcast-compatible with *areal_density*) of pixel
        areas in cm².
    mask:
        Boolean array of the same shape as *areal_density*; True pixels are
        included in the sum.

    Returns
    -------
    float
        Total mass in grams.

    Raises
    ------
    ValueError
        If *pixel_area_cm2* is invalid (None, zero, or negative), or if *mask*
        contains no True pixels.
    """
    validate_pixel_area(pixel_area_cm2)
    if not np.any(mask):
        raise ValueError(
            "mask has no True pixels; cannot compute regional mass for an empty region. "
            "Ensure the mask was generated correctly and spatially corresponds to this image."
        )
    area = np.broadcast_to(np.asarray(pixel_area_cm2, dtype=np.float64), areal_density.shape)
    return float(np.nansum(areal_density[mask].astype(np.float64) * area[mask]))


def center_of_mass_mm(
    areal_density: np.ndarray,
    pixel_area_cm2: np.ndarray | float,
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    mask: np.ndarray,
) -> tuple[float, float]:
    """Mass-weighted centroid in mm.

    Parameters
    ----------
    areal_density:
        2-D array of areal density values in g/cm².
    pixel_area_cm2:
        Scalar or array of pixel areas in cm².
    x_mm, y_mm:
        Coordinate arrays of shape (rows, cols) in mm.
        x is lateral, y is superior-inferior.
    mask:
        Boolean inclusion mask.

    Returns
    -------
    tuple[float, float]
        ``(x_com_mm, y_com_mm)`` in mm.

    Raises
    ------
    ValueError
        If total mass in the masked region is zero or pixel area is invalid.
    """
    validate_pixel_area(pixel_area_cm2)
    area = np.broadcast_to(np.asarray(pixel_area_cm2, dtype=np.float64), areal_density.shape)
    dm = areal_density.astype(np.float64) * area
    total_mass = float(np.nansum(dm[mask]))
    if total_mass == 0.0:
        raise ValueError(
            "Total mass in the masked region is zero; center of mass is undefined. "
            "Check that areal_density contains non-zero values within the masked region."
        )
    x_com = float(np.nansum(x_mm[mask] * dm[mask])) / total_mass
    y_com = float(np.nansum(y_mm[mask] * dm[mask])) / total_mass
    return x_com, y_com


def planar_moment_of_inertia(
    areal_density: np.ndarray,
    pixel_area_cm2: np.ndarray | float,
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    mask: np.ndarray,
    cx_mm: float | None = None,
    cy_mm: float | None = None,
) -> float:
    """Moment of inertia about the axis normal to the detector (g·mm²).

    Parameters
    ----------
    areal_density:
        2-D array of areal density values in g/cm².
    pixel_area_cm2:
        Scalar or array of pixel areas in cm².
    x_mm, y_mm:
        Coordinate arrays of shape (rows, cols) in mm.
    mask:
        Boolean inclusion mask.
    cx_mm, cy_mm:
        Pivot point in mm.  If either is None, the center of mass is computed
        and used as the pivot.

    Returns
    -------
    float
        Planar moment of inertia in g·mm².
    """
    validate_pixel_area(pixel_area_cm2)
    if cx_mm is None or cy_mm is None:
        cx_mm, cy_mm = center_of_mass_mm(
            areal_density, pixel_area_cm2, x_mm, y_mm, mask
        )
    area = np.broadcast_to(np.asarray(pixel_area_cm2, dtype=np.float64), areal_density.shape)
    dm = areal_density.astype(np.float64) * area
    r_sq = (x_mm - cx_mm) ** 2 + (y_mm - cy_mm) ** 2
    return float(np.nansum(dm[mask] * r_sq[mask]))


def compute_mass_properties(
    bundle: "QuantitativeMapBundle",
    region_masks: dict[str, np.ndarray] | None = None,
) -> dict[str, dict]:
    """Compute mass, CoM, and I_perp for whole body and each provided region mask.

    Parameters
    ----------
    bundle:
        ``QuantitativeMapBundle`` containing maps, masks, and geometry.
    region_masks:
        Optional additional named boolean masks keyed by region name.
        Merged with the bundle's whole-body valid mask under the key
        ``"whole_body"``.

    Returns
    -------
    dict[str, dict]
        Nested dict::

            {
                region_name: {
                    "fat_g": float,
                    "lean_g": float,
                    "bmc_g": float,
                    "total_g": float,
                    "x_cm_mm": float,
                    "y_cm_mm": float,
                    "I_perp_g_mm2": float,
                }
            }

        Components that are unavailable in the bundle are stored as ``nan``.

    Raises
    ------
    ValueError
        If pixel_area is unavailable in the bundle geometry, or if no density
        map is present for mechanical property computation.
    """
    validate_pixel_area(bundle.geometry.pixel_area_cm2)

    pixel_area = bundle.geometry.pixel_area_cm2
    x_mm: np.ndarray = bundle.geometry.x_coordinates_mm
    y_mm: np.ndarray = bundle.geometry.y_coordinates_mm

    valid_mask: np.ndarray = (
        bundle.masks.valid
        if bundle.masks.valid is not None
        else bundle.masks.body
    )

    all_masks: dict[str, np.ndarray] = {"whole_body": valid_mask}
    if region_masks:
        all_masks.update(region_masks)

    fat: np.ndarray | None = bundle.maps.fat_areal_density_g_cm2
    lean: np.ndarray | None = bundle.maps.lean_areal_density_g_cm2
    bmc: np.ndarray | None = bundle.maps.bmc_areal_density_g_cm2
    total: np.ndarray | None = bundle.maps.total_areal_density_g_cm2

    density_for_mechanics: np.ndarray | None = next(
        (arr for arr in (total, fat, lean, bmc) if arr is not None), None
    )
    if density_for_mechanics is None:
        raise ValueError(
            "No density map is available in the bundle (fat, lean, bmc, and total are all "
            "None). Load at least one quantitative map before computing mass properties."
        )

    results: dict[str, dict] = {}

    for region_name, mask in all_masks.items():
        if not np.any(mask):
            log.warning("empty_mask_skipped", region=region_name)
            continue

        fat_g = regional_mass_g(fat, pixel_area, mask) if fat is not None else float("nan")
        lean_g = (
            regional_mass_g(lean, pixel_area, mask) if lean is not None else float("nan")
        )
        bmc_g = regional_mass_g(bmc, pixel_area, mask) if bmc is not None else float("nan")

        if total is not None:
            total_g = regional_mass_g(total, pixel_area, mask)
        else:
            components = [v for v in (fat_g, lean_g, bmc_g) if not np.isnan(v)]
            total_g = float(sum(components)) if components else float("nan")

        x_com, y_com = center_of_mass_mm(
            density_for_mechanics, pixel_area, x_mm, y_mm, mask
        )
        i_perp = planar_moment_of_inertia(
            density_for_mechanics, pixel_area, x_mm, y_mm, mask, x_com, y_com
        )

        results[region_name] = {
            "fat_g": fat_g,
            "lean_g": lean_g,
            "bmc_g": bmc_g,
            "total_g": total_g,
            "x_cm_mm": x_com,
            "y_cm_mm": y_com,
            "I_perp_g_mm2": i_perp,
        }
        log.info(
            "mass_properties_computed",
            region=region_name,
            total_g=round(total_g, 1) if not np.isnan(total_g) else None,
            x_cm_mm=round(x_com, 2),
            y_cm_mm=round(y_com, 2),
        )

    return results
