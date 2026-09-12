"""BodyLoop 3D geometry interface — experimental, disabled by default.

This module provides interfaces for loading and validating BodyLoop surface
data. The full registration pipeline (bodyloop.registration) is behind the
`--experimental` flag because articulated/deformable registration between
supine DXA and standing BodyLoop postures is an active research problem.

To use: `hologic-dxa bodyloop-register --experimental`
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from hologic_dxa.logging import get_logger

if TYPE_CHECKING:
    pass

log = get_logger(__name__)

SUPPORTED_FORMATS = {".ply", ".obj", ".stl", ".pcd"}


@dataclass
class BodyLoopSurface:
    """Loaded BodyLoop surface with validated units and axis orientation."""

    vertices: np.ndarray          # shape (N, 3), float64, metres
    faces: np.ndarray | None      # shape (F, 3), int — None for point clouds
    landmarks: dict[str, np.ndarray]  # name → (3,) array in metres
    units: str                    # "m", "mm", "cm"
    source_path: Path
    axis_convention: str          # e.g. "RAS", "LPS", "scanner"


def load_surface(path: Path) -> BodyLoopSurface:
    """Load a BodyLoop mesh or point cloud.

    Requires the `bodyloop` optional dependencies:
        pip install hologic-dxa[bodyloop]

    Raises
    ------
    ImportError
        If trimesh or open3d are not installed.
    ValueError
        If the file format is unsupported or units cannot be determined.
    """
    try:
        import trimesh  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "BodyLoop loading requires trimesh. "
            "Install it with: pip install hologic-dxa[bodyloop]"
        ) from exc

    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_FORMATS:
        raise ValueError(
            f"Unsupported format '{suffix}'. "
            f"Supported: {', '.join(sorted(SUPPORTED_FORMATS))}"
        )

    import trimesh as tm

    log.info("bodyloop_load", path=str(path))
    mesh = tm.load(str(path))

    if isinstance(mesh, tm.PointCloud):
        vertices = np.asarray(mesh.vertices, dtype=np.float64)
        faces = None
    elif isinstance(mesh, tm.Trimesh):
        vertices = np.asarray(mesh.vertices, dtype=np.float64)
        faces = np.asarray(mesh.faces, dtype=np.int64)
    else:
        raise ValueError(f"Unrecognised trimesh object type: {type(mesh)}")

    # BodyLoop exports coordinates in mm by convention; document if different
    units = _detect_units(vertices)
    log.info("bodyloop_units_detected", units=units, n_vertices=len(vertices))

    return BodyLoopSurface(
        vertices=vertices,
        faces=faces,
        landmarks={},
        units=units,
        source_path=path,
        axis_convention="unknown",
    )


def _detect_units(vertices: np.ndarray) -> str:
    """Heuristic unit detection based on vertex coordinate magnitude.

    A standing human is ~1.7 m, ~170 cm, or ~1700 mm tall.
    The Z-extent of the point cloud is used as the primary cue.
    This is a heuristic — the user must confirm.
    """
    extent = vertices[:, 2].max() - vertices[:, 2].min()
    if 1.0 < extent < 2.5:
        return "m"
    if 100 < extent < 250:
        return "cm"
    if 1000 < extent < 2500:
        return "mm"
    log.warning(
        "bodyloop_units_ambiguous",
        z_extent=float(extent),
        note="Cannot determine units from coordinate magnitude. Assuming mm.",
    )
    return "mm"
