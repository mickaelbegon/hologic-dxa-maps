"""BodyLoop ↔ DXA registration — experimental, research-grade.

The registration problem between a standing BodyLoop surface and a supine
DXA scan is fundamentally non-rigid:

1. Gravity deforms soft tissue differently in the two postures.
2. The table creates dorsal flattening in the DXA scan.
3. Limb positions differ (arms alongside the body in DXA).
4. Spine curvature changes between standing and supine.

A rigid or even articulated-rigid registration is therefore only a first
approximation. This module is provided as an experimental interface to be
refined as data and algorithms become available.

This module is DISABLED by default. Activate with:
    hologic-dxa bodyloop-register --experimental
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from hologic_dxa.bodyloop.interface import BodyLoopSurface


@dataclass
class RegistrationResult:
    """Result of a BodyLoop ↔ DXA registration attempt."""

    transform_4x4: np.ndarray    # homogeneous 4×4 rigid transform (initial estimate)
    residual_mm: float           # RMS landmark residual after registration
    method: str                  # e.g. "landmark_rigid", "icp_articulated"
    converged: bool
    n_iterations: int
    warnings: list[str]


def rigid_landmark_registration(
    bodyloop_landmarks: dict[str, np.ndarray],
    dxa_landmarks: dict[str, np.ndarray],
) -> RegistrationResult:
    """Estimate a rigid transform from anatomical landmark correspondences.

    Parameters
    ----------
    bodyloop_landmarks:
        Landmark coordinates in BodyLoop space (mm).
    dxa_landmarks:
        Corresponding landmark coordinates in DXA image space (mm).

    Returns
    -------
    RegistrationResult
        The transform is an INITIAL ESTIMATE ONLY. Do not report it as a
        validated registration without deformable refinement.

    Raises
    ------
    NotImplementedError
        Always — this function is a placeholder until landmark data and
        a validated algorithm are available.
    """
    raise NotImplementedError(
        "Landmark-based rigid registration is not yet implemented. "
        "Required inputs: (1) BodyLoop anatomical landmarks in mm, "
        "(2) corresponding DXA image landmarks in mm, "
        "(3) a validated algorithm accounting for posture change. "
        "See docs/limitations.md#bodyloop-registration."
    )


def deformable_registration(
    surface: "BodyLoopSurface",
    dxa_silhouette: np.ndarray,
    initial_transform: np.ndarray | None = None,
) -> RegistrationResult:
    """Deformable registration of BodyLoop surface to DXA silhouette.

    This requires a non-rigid deformation model that accounts for:
    - Gravity-induced tissue redistribution (supine vs. standing)
    - Table contact deformation (dorsal flattening)
    - Limb repositioning

    Raises
    ------
    NotImplementedError
        Always — this is an open research problem.
    """
    raise NotImplementedError(
        "Deformable BodyLoop ↔ DXA registration is an open research problem. "
        "No validated algorithm is implemented. "
        "Contributions are welcome — see CONTRIBUTING.md."
    )
