"""Fan-beam ray projection utilities for BodyLoop thickness maps.

Given a BodyLoop mesh articulated to match the DXA posture, this module
computes the line integral of the body indicator function along each
DXA detector ray:

    T_{ij} = ∫_{R_{ij}} 1_body(r) ds

This thickness map can then be combined with the DXA areal densities to
estimate the volumetric density distribution.

Status: EXPERIMENTAL — requires a registered mesh (bodyloop.registration)
and validated DXA source geometry (source-to-detector distance, pixel size).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from hologic_dxa.bodyloop.interface import BodyLoopSurface
    from hologic_dxa.maps.bundle import QuantitativeMapBundle


@dataclass
class ThicknessMap:
    """Pixel-level thickness of the body along DXA projection rays (mm)."""

    values_mm: np.ndarray          # shape (rows, cols), float64
    valid_mask: np.ndarray         # bool, True where ray intersected body
    source_position_mm: np.ndarray  # (3,) — DXA X-ray source
    detector_normal_mm: np.ndarray  # (3,) — normal to detector plane


def compute_thickness_map(
    surface: "BodyLoopSurface",
    bundle: "QuantitativeMapBundle",
    source_position_mm: np.ndarray,
) -> ThicknessMap:
    """Compute ray-intersection thickness map from a registered BodyLoop mesh.

    Raises
    ------
    NotImplementedError
        Until a validated registered mesh and DXA geometry are available.
    ImportError
        If trimesh is not installed.
    """
    try:
        import trimesh  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Thickness map computation requires trimesh. "
            "Install with: pip install hologic-dxa[bodyloop]"
        ) from exc

    raise NotImplementedError(
        "Thickness map computation requires: "
        "(1) a validated registered BodyLoop mesh (see bodyloop.registration), "
        "(2) confirmed DXA source geometry (source position, SDD), "
        "(3) a closed watertight mesh for ray intersection. "
        "This is currently blocked on registration validation."
    )
