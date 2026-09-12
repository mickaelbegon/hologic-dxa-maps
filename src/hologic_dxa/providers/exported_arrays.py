"""ExportedArraysProvider — loads pre-calibrated tissue maps from a local directory.

Expected layout
---------------
fat.npy | fat.csv
lean.npy | lean.csv
bmc.npy | bmc.csv
valid_mask.npy          (optional; defaults to all-True)
metadata.json

metadata.json schema (required fields)
---------------------------------------
{
    "units":                     {"fat": "g/cm2", "lean": "g/cm2", "bmc": "g/cm2"},
    "pixel_area_cm2":            <float > 0>,
    "source_sop_instance_uid":   "<non-empty string>",
    "already_geometry_corrected": <bool>
}
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import structlog

from hologic_dxa.maps.bundle import QuantitativeMapBundle
from hologic_dxa.providers.base import QuantitativeDataUnavailableError

log = structlog.get_logger(__name__)

_REQUIRED_UNIT = "g/cm2"
_COMPONENTS = ("fat", "lean", "bmc")


class ExportedArraysProvider:
    name = "exported_arrays"
    experimental = False

    def can_read(self, source: Path) -> bool:
        if not source.is_dir():
            return False
        has_metadata = (source / "metadata.json").is_file()
        has_any_component = any(
            (source / f"{comp}.npy").is_file() or (source / f"{comp}.csv").is_file()
            for comp in _COMPONENTS
        )
        return has_metadata and has_any_component

    def load(self, source: Path) -> QuantitativeMapBundle:
        log.info("exported_arrays_load_start")

        metadata = _load_metadata(source)
        arrays = {comp: _load_array(source, comp) for comp in _COMPONENTS}
        _validate_units(metadata)
        _validate_pixel_area(metadata)
        _validate_sop_uid(metadata)
        _validate_shapes(arrays)
        valid_mask = _load_valid_mask(source, reference_shape=arrays["fat"].shape)
        _validate_values(arrays, valid_mask)

        log.info(
            "exported_arrays_load_ok",
            shape=arrays["fat"].shape,
            pixel_area_cm2=metadata["pixel_area_cm2"],
        )
        return QuantitativeMapBundle(
            fat=arrays["fat"],
            lean=arrays["lean"],
            bmc=arrays["bmc"],
            valid_mask=valid_mask,
            pixel_area_cm2=float(metadata["pixel_area_cm2"]),
            source_sop_instance_uid=str(metadata["source_sop_instance_uid"]),
            already_geometry_corrected=bool(metadata.get("already_geometry_corrected", True)),
            provider_name=self.name,
            metadata=metadata,
        )

    def describe_capabilities(self) -> dict[str, object]:
        return {
            "status": "available",
            "description": (
                "Reads pre-calibrated fat/lean/bmc arrays (npy or csv) plus "
                "metadata.json from a local directory.  Units must be g/cm2.  "
                "total = fat + lean + bmc is always derived, never stored."
            ),
            "required_files": [
                "fat.npy or fat.csv",
                "lean.npy or lean.csv",
                "bmc.npy or bmc.csv",
                "metadata.json",
            ],
            "optional_files": ["valid_mask.npy"],
            "accepted_units": _REQUIRED_UNIT,
            "experimental": False,
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_metadata(source: Path) -> dict[str, object]:
    meta_path = source / "metadata.json"
    if not meta_path.is_file():
        raise QuantitativeDataUnavailableError(
            f"metadata.json not found in '{source}'. "
            "This file is required and must contain: units, pixel_area_cm2, "
            "source_sop_instance_uid, already_geometry_corrected."
        )
    with meta_path.open(encoding="utf-8") as fh:
        return json.load(fh)  # type: ignore[no-any-return]


def _load_array(source: Path, component: str) -> np.ndarray:
    npy = source / f"{component}.npy"
    csv = source / f"{component}.csv"

    if npy.is_file():
        return np.load(npy)

    if csv.is_file():
        return np.loadtxt(csv, delimiter=",")

    raise QuantitativeDataUnavailableError(
        f"Component '{component}' not found in '{source}'. "
        f"Provide either '{component}.npy' or '{component}.csv'."
    )


def _load_valid_mask(source: Path, reference_shape: tuple[int, ...]) -> np.ndarray:
    mask_path = source / "valid_mask.npy"
    if mask_path.is_file():
        mask = np.load(mask_path).astype(bool)
        if mask.shape != reference_shape:
            raise QuantitativeDataUnavailableError(
                f"valid_mask.npy shape {mask.shape} does not match "
                f"component shape {reference_shape}."
            )
        return mask
    return np.ones(reference_shape, dtype=bool)


def _validate_units(metadata: dict[str, object]) -> None:
    units = metadata.get("units")
    if not isinstance(units, dict):
        raise QuantitativeDataUnavailableError(
            "metadata.json is missing the 'units' dict. "
            f"Expected {{\"fat\": \"{_REQUIRED_UNIT}\", ...}}."
        )
    bad = {
        comp: units.get(comp)
        for comp in _COMPONENTS
        if units.get(comp) != _REQUIRED_UNIT
    }
    if bad:
        raise QuantitativeDataUnavailableError(
            f"Unit mismatch in metadata.json: {bad}. "
            f"Only '{_REQUIRED_UNIT}' is accepted without explicit conversion. "
            "Re-export the arrays with the correct units."
        )


def _validate_pixel_area(metadata: dict[str, object]) -> None:
    area = metadata.get("pixel_area_cm2")
    try:
        area_f = float(area)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise QuantitativeDataUnavailableError(
            "metadata.json 'pixel_area_cm2' must be a positive float."
        ) from exc
    if area_f <= 0:
        raise QuantitativeDataUnavailableError(
            f"pixel_area_cm2 must be > 0, got {area_f}. "
            "Check the scanner calibration or metadata export."
        )


def _validate_sop_uid(metadata: dict[str, object]) -> None:
    uid = metadata.get("source_sop_instance_uid")
    if not uid or not str(uid).strip():
        raise QuantitativeDataUnavailableError(
            "metadata.json 'source_sop_instance_uid' is absent or empty. "
            "This field is required for provenance tracking."
        )


def _validate_shapes(arrays: dict[str, np.ndarray]) -> None:
    shapes = {comp: arr.shape for comp, arr in arrays.items()}
    unique = set(shapes.values())
    if len(unique) != 1:
        raise QuantitativeDataUnavailableError(
            f"All component arrays must share the same shape, got: {shapes}. "
            "Re-export so fat, lean, and bmc arrays have identical dimensions."
        )


def _validate_values(arrays: dict[str, np.ndarray], valid_mask: np.ndarray) -> None:
    _FLOOR = -1e-6
    for comp, arr in arrays.items():
        masked_vals = arr[valid_mask]
        if masked_vals.size == 0:
            continue
        min_val = float(masked_vals.min())
        if min_val < _FLOOR:
            raise QuantitativeDataUnavailableError(
                f"Component '{comp}' contains values as low as {min_val:.6g} g/cm2 "
                "in the valid region (threshold: -1e-6). "
                "Negative tissue density indicates a calibration or export error."
            )
