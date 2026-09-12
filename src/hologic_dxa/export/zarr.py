"""Zarr export for quantitative map bundles.

Zarr is an alternative to HDF5 with better cloud storage support
(chunked, compressed arrays stored as directories or zip archives).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from hologic_dxa.logging import get_logger

if TYPE_CHECKING:
    from hologic_dxa.models import QuantitativeMapBundle

log = get_logger(__name__)


def export_bundle_zarr(
    bundle: "QuantitativeMapBundle",
    output_path: Path,
    compressor: object | None = None,
) -> None:
    """Export QuantitativeMapBundle to a Zarr store.

    Parameters
    ----------
    output_path:
        Path to the Zarr store directory (created if absent).
    compressor:
        numcodecs compressor instance. None uses Zarr's default (Blosc).

    Raises
    ------
    ImportError
        If zarr is not installed.
    """
    try:
        import zarr  # type: ignore[import-untyped]
    except ImportError as exc:
        raise ImportError(
            "Zarr export requires zarr. Install with: pip install zarr"
        ) from exc

    import numpy as np

    store = zarr.open(str(output_path), mode="w")
    maps_group = store.require_group("maps")
    masks_group = store.require_group("masks")
    geom_group = store.require_group("geometry")
    meta_group = store.require_group("metadata")

    _write_map(maps_group, "fat_areal_density_g_cm2", bundle.fat, compressor)
    _write_map(maps_group, "lean_areal_density_g_cm2", bundle.lean, compressor)
    _write_map(maps_group, "bmc_areal_density_g_cm2", bundle.bmc, compressor)
    _write_map(maps_group, "total_areal_density_g_cm2", bundle.total, compressor)

    if bundle.body_mask is not None:
        masks_group.array("body", bundle.body_mask, dtype=bool, compressor=compressor)
    if bundle.bone_mask is not None:
        masks_group.array("bone", bundle.bone_mask, dtype=bool, compressor=compressor)
    masks_group.array("valid", bundle.valid_mask, dtype=bool, compressor=compressor)

    geom = bundle.geometry
    if isinstance(geom.pixel_spacing_mm, (list, tuple)):
        geom_group.attrs["pixel_spacing_mm"] = list(geom.pixel_spacing_mm)
    geom_group.attrs["rows"] = geom.rows
    geom_group.attrs["columns"] = geom.columns
    geom_group.attrs["is_fan_beam_corrected"] = geom.is_fan_beam_corrected

    meta_group.attrs.update(bundle.metadata)
    meta_group.attrs["pipeline_version"] = bundle.provenance.pipeline_version
    meta_group.attrs["timestamp_utc"] = bundle.provenance.timestamp_utc.isoformat()

    log.info("zarr_export_complete", path=str(output_path))


def _write_map(
    group: object,
    name: str,
    qmap: object | None,
    compressor: object | None,
) -> None:
    import zarr  # type: ignore[import-untyped]

    if qmap is None:
        return
    arr = getattr(group, "array", None)
    if arr is None:
        return
    dataset = arr(name, qmap.values, dtype="float64", compressor=compressor)
    dataset.attrs["unit"] = qmap.unit
    dataset.attrs["quantity"] = qmap.quantity
