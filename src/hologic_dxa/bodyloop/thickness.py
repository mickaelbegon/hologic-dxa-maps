"""Thickness map post-processing utilities."""

from __future__ import annotations

import numpy as np

from hologic_dxa.bodyloop.projection import ThicknessMap


def smooth_thickness_map(
    thickness: ThicknessMap,
    sigma_mm: float = 5.0,
    pixel_size_mm: float = 1.0,
) -> ThicknessMap:
    """Apply Gaussian smoothing to the thickness map.

    Parameters
    ----------
    sigma_mm:
        Standard deviation of the Gaussian kernel in mm.
    pixel_size_mm:
        Pixel size in mm (used to convert sigma_mm to pixels).
    """
    from scipy.ndimage import gaussian_filter  # type: ignore[import-untyped]

    sigma_px = sigma_mm / pixel_size_mm
    smoothed = gaussian_filter(thickness.values_mm, sigma=sigma_px)
    smoothed[~thickness.valid_mask] = 0.0

    return ThicknessMap(
        values_mm=smoothed,
        valid_mask=thickness.valid_mask,
        source_position_mm=thickness.source_position_mm,
        detector_normal_mm=thickness.detector_normal_mm,
    )


def volumetric_density_estimate(
    areal_density_g_cm2: np.ndarray,
    thickness_mm: np.ndarray,
    valid_mask: np.ndarray,
) -> np.ndarray:
    """Estimate volumetric density from areal density and thickness.

    volumetric_density [g/cm³] = areal_density [g/cm²] / thickness [cm]

    This estimate is only valid where:
    - The thickness is non-zero
    - The tissue is approximately uniform along the ray
    - Fan-beam geometry corrections have been applied

    Returns nan where thickness ≤ 0 or where valid_mask is False.
    """
    thickness_cm = thickness_mm / 10.0
    result = np.full_like(areal_density_g_cm2, np.nan)
    valid = valid_mask & (thickness_cm > 0)
    result[valid] = areal_density_g_cm2[valid] / thickness_cm[valid]
    return result
