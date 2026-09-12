"""HDF5 export for QuantitativeMapBundle.

Produces a self-describing HDF5 file with units attributes on every dataset,
provenance tracking, and graceful handling of missing map components.

HDF5 structure
--------------
/maps/fat_areal_density_g_cm2         (float32, units=g/cm^2)
/maps/lean_areal_density_g_cm2        (float32, units=g/cm^2)
/maps/bmc_areal_density_g_cm2         (float32, units=g/cm^2)
/maps/total_areal_density_g_cm2       (float32, units=g/cm^2)
/masks/body                            (uint8/bool)
/masks/bone                            (uint8/bool)
/masks/valid                           (uint8/bool)
/geometry/pixel_area_cm2              (float64, units=cm^2)
/geometry/x_coordinates_mm            (float64, units=mm)
/geometry/y_coordinates_mm            (float64, units=mm)
/metadata/                             (group; key metadata stored as attrs)
/provenance/input_hashes              (JSON string dataset)
/provenance/operations                (JSON string dataset)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np

from hologic_dxa import __version__
from hologic_dxa.logging import get_logger

if TYPE_CHECKING:
    from hologic_dxa.models import QuantitativeMapBundle

log = get_logger(__name__)


def export_bundle_hdf5(
    bundle: "QuantitativeMapBundle",
    output_path: Path,
    compression: str = "gzip",
    compression_opts: int = 4,
) -> None:
    """Export QuantitativeMapBundle to HDF5.

    Parameters
    ----------
    bundle:
        Bundle to export.
    output_path:
        Destination ``.h5`` file.  Parent directories are created if absent.
    compression:
        HDF5 compression filter (default ``"gzip"``).
    compression_opts:
        Compression level (default 4; gzip range 0–9).

    Notes
    -----
    Missing datasets are skipped but listed in the
    ``/metadata`` group attribute ``missing_datasets`` (JSON array).
    All datasets carry a ``units`` HDF5 attribute.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    compress_kwargs: dict = {
        "compression": compression,
        "compression_opts": compression_opts,
    }

    missing: list[str] = []

    # --- Resolve geometry-derived arrays -----------------------------------
    # bundle.geometry is a MapGeometry (rows, columns, pixel_spacing_mm, origin_mm).
    # Coordinate arrays and pixel area may be stored on individual QuantitativeMap
    # objects; fall back to computing them from MapGeometry.
    from hologic_dxa.maps.geometry import pixel_coordinates_mm as _pixel_coords_mm

    first_qmap = next(
        (m for m in (bundle.fat, bundle.lean, bundle.bmc, bundle.total) if m is not None),
        None,
    )

    if first_qmap is not None and first_qmap.pixel_area_cm2 is not None:
        pixel_area_data: np.ndarray | float = first_qmap.pixel_area_cm2
    else:
        pixel_area_data = bundle.geometry.pixel_area_cm2_scalar()

    if first_qmap is not None and first_qmap.x_coordinates_mm is not None:
        x_coords: np.ndarray | None = first_qmap.x_coordinates_mm
        y_coords: np.ndarray | None = first_qmap.y_coordinates_mm
    else:
        x_coords, y_coords = _pixel_coords_mm(
            bundle.geometry.rows,
            bundle.geometry.columns,
            bundle.geometry.pixel_spacing_mm[0],
            bundle.geometry.pixel_spacing_mm[1],
            bundle.geometry.origin_mm,
        )

    with h5py.File(output_path, "w") as f:
        # ------------------------------------------------------------------ maps
        maps_grp = f.create_group("maps")
        _write_float_dataset(
            maps_grp, "fat_areal_density_g_cm2",
            bundle.fat.values if bundle.fat is not None else None,
            units="g/cm^2", missing=missing, **compress_kwargs,
        )
        _write_float_dataset(
            maps_grp, "lean_areal_density_g_cm2",
            bundle.lean.values if bundle.lean is not None else None,
            units="g/cm^2", missing=missing, **compress_kwargs,
        )
        _write_float_dataset(
            maps_grp, "bmc_areal_density_g_cm2",
            bundle.bmc.values if bundle.bmc is not None else None,
            units="g/cm^2", missing=missing, **compress_kwargs,
        )
        _write_float_dataset(
            maps_grp, "total_areal_density_g_cm2",
            bundle.total.values if bundle.total is not None else None,
            units="g/cm^2", missing=missing, **compress_kwargs,
        )

        # ----------------------------------------------------------------- masks
        masks_grp = f.create_group("masks")
        _write_bool_dataset(
            masks_grp, "body",
            getattr(bundle, "body_mask", None),
            missing=missing, **compress_kwargs,
        )
        _write_bool_dataset(
            masks_grp, "bone",
            getattr(bundle, "bone_mask", None),
            missing=missing, **compress_kwargs,
        )
        _write_bool_dataset(
            masks_grp, "valid",
            bundle.valid_mask,
            missing=missing, **compress_kwargs,
        )

        # --------------------------------------------------------------- geometry
        geom_grp = f.create_group("geometry")
        _write_float_dataset(
            geom_grp, "pixel_area_cm2",
            pixel_area_data,
            units="cm^2", missing=missing, dtype=np.float64, **compress_kwargs,
        )
        _write_float_dataset(
            geom_grp, "x_coordinates_mm",
            x_coords,
            units="mm", missing=missing, dtype=np.float64, **compress_kwargs,
        )
        _write_float_dataset(
            geom_grp, "y_coordinates_mm",
            y_coords,
            units="mm", missing=missing, dtype=np.float64, **compress_kwargs,
        )

        # -------------------------------------------------------------- metadata
        meta_grp = f.create_group("metadata")
        prov = getattr(bundle, "provenance", None)
        _write_metadata_attrs(meta_grp, bundle, prov)
        if missing:
            meta_grp.attrs["missing_datasets"] = json.dumps(missing)

        # ------------------------------------------------------------- provenance
        prov_grp = f.create_group("provenance")
        if prov is not None:
            input_hashes = getattr(prov, "input_hashes", None) or {}
            operations = getattr(prov, "operations", None) or []
            prov_grp.create_dataset(
                "input_hashes",
                data=json.dumps(input_hashes),
            )
            prov_grp.create_dataset(
                "operations",
                data=json.dumps(operations),
            )
        else:
            missing.extend(["provenance/input_hashes", "provenance/operations"])
            log.warning("hdf5_provenance_missing", path=str(output_path))

        # Backfill missing list in case provenance was appended after metadata group
        if missing and "missing_datasets" not in meta_grp.attrs:
            meta_grp.attrs["missing_datasets"] = json.dumps(missing)

    log.info(
        "hdf5_exported",
        path=str(output_path),
        missing_datasets=missing,
        compression=compression,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _write_float_dataset(
    group: h5py.Group,
    name: str,
    data: np.ndarray | float | None,
    *,
    units: str,
    missing: list[str],
    dtype: type = np.float32,
    compression: str = "gzip",
    compression_opts: int = 4,
) -> None:
    if data is None:
        missing.append(f"{group.name}/{name}")
        log.debug("hdf5_dataset_skipped", group=group.name, name=name)
        return
    arr = np.asarray(data, dtype=dtype)
    # Scalar pixel_area_cm2 stored as 0-D dataset; disable chunked compression for scalars
    if arr.ndim == 0:
        ds = group.create_dataset(name, data=arr)
    else:
        ds = group.create_dataset(
            name, data=arr, compression=compression, compression_opts=compression_opts
        )
    ds.attrs["units"] = units


def _write_bool_dataset(
    group: h5py.Group,
    name: str,
    data: np.ndarray | None,
    *,
    missing: list[str],
    compression: str = "gzip",
    compression_opts: int = 4,
) -> None:
    if data is None:
        missing.append(f"{group.name}/{name}")
        log.debug("hdf5_dataset_skipped", group=group.name, name=name)
        return
    arr = np.asarray(data, dtype=np.uint8)
    ds = group.create_dataset(
        name, data=arr, compression=compression, compression_opts=compression_opts
    )
    ds.attrs["units"] = "boolean (0=False, 1=True)"


def _write_metadata_attrs(
    group: h5py.Group,
    bundle: "QuantitativeMapBundle",
    prov: object | None,
) -> None:
    """Write device and pipeline metadata as HDF5 group attributes."""
    group.attrs["pipeline_version"] = __version__

    if prov is not None:
        for attr in (
            "device_manufacturer",
            "device_model",
            "apex_version",
            "source_sop_instance_uid",
            "study_instance_uid",
            "series_instance_uid",
            "transfer_syntax_uid",
        ):
            val = getattr(prov, attr, None)
            if val is not None:
                group.attrs[attr] = str(val)

        ts = getattr(prov, "timestamp_utc", None)
        if ts is not None:
            group.attrs["provenance_timestamp_utc"] = ts.isoformat()

    # Merge any extra metadata the bundle exposes directly
    extra_meta = getattr(bundle, "metadata", None)
    if isinstance(extra_meta, dict):
        for k, v in extra_meta.items():
            if k not in group.attrs:
                try:
                    group.attrs[k] = str(v)
                except Exception:  # noqa: BLE001
                    pass
