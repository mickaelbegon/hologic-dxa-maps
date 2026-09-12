"""Domain model classes for DXA quantitative maps.

All array-holding classes use ``dataclasses.dataclass`` because pydantic v2
cannot validate numpy arrays. Pydantic is reserved for config/settings models
(see :mod:`hologic_dxa.config`).

Physical notes
--------------
Fan-beam DXA scanners produce images where the pixel area varies across the
field of view due to beam divergence. Scalar pixel-area approximations are
only valid for flat-detector (pencil-beam) geometry. Always prefer the
per-pixel ``pixel_area_cm2`` array when it is available.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from hologic_dxa.provenance import Provenance


# ---------------------------------------------------------------------------
# QuantitativeMap
# ---------------------------------------------------------------------------

@dataclass
class QuantitativeMap:
    """A single physical quantity sampled on a 2-D pixel grid.

    Parameters
    ----------
    values:
        2-D float64 array of physical quantities (shape ``(rows, cols)``).
    quantity:
        Machine-readable quantity name, e.g. ``"fat_areal_density"``.
    unit:
        Unit string compatible with pint, e.g. ``"g/cm2"``.
    valid_mask:
        Boolean array (same shape as *values*). ``True`` where the pixel is
        anatomically valid and the value should be trusted.
    pixel_area_cm2:
        Scalar or per-pixel area in cm². ``None`` if unknown.
    x_coordinates_mm:
        Column-wise x-positions of pixel centres (mm). May be ``None``.
    y_coordinates_mm:
        Row-wise y-positions of pixel centres (mm). May be ``None``.
    source_uid:
        SOP Instance UID of the originating DICOM object.
    provenance:
        JSON-serialisable provenance dict (from :meth:`Provenance.to_dict`).
    uncertainty:
        Optional per-pixel uncertainty array (same shape as *values*).
    """

    values: np.ndarray
    quantity: str
    unit: str
    valid_mask: np.ndarray
    pixel_area_cm2: np.ndarray | float | None
    x_coordinates_mm: np.ndarray | None
    y_coordinates_mm: np.ndarray | None
    source_uid: str
    provenance: dict[str, Any]
    uncertainty: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.values.ndim != 2:
            raise ValueError(
                f"QuantitativeMap.values must be a 2-D array, got shape {self.values.shape}."
            )
        if self.values.dtype != np.float64:
            self.values = self.values.astype(np.float64)

        if self.valid_mask.shape != self.values.shape:
            raise ValueError(
                f"valid_mask shape {self.valid_mask.shape} does not match "
                f"values shape {self.values.shape}."
            )
        if self.valid_mask.dtype != bool:
            self.valid_mask = self.valid_mask.astype(bool)

        if isinstance(self.pixel_area_cm2, np.ndarray):
            if self.pixel_area_cm2.shape != self.values.shape:
                raise ValueError(
                    f"pixel_area_cm2 array shape {self.pixel_area_cm2.shape} does not "
                    f"match values shape {self.values.shape}."
                )

        if self.uncertainty is not None and self.uncertainty.shape != self.values.shape:
            raise ValueError(
                f"uncertainty shape {self.uncertainty.shape} does not match "
                f"values shape {self.values.shape}."
            )

        if not self.quantity:
            raise ValueError("quantity must be a non-empty string.")
        if not self.unit:
            raise ValueError("unit must be a non-empty string.")
        if not self.source_uid:
            raise ValueError("source_uid must be a non-empty string.")

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------

    def integrated_mass_g(self) -> float:
        """Sum of values × pixel_area over all valid pixels, in grams.

        Raises
        ------
        ValueError
            If ``pixel_area_cm2`` is ``None`` (mass calculation requires area).
        """
        if self.pixel_area_cm2 is None:
            raise ValueError(
                "Cannot compute integrated mass: pixel_area_cm2 is None. "
                "Supply a scalar or per-pixel area array at construction time."
            )
        masked_values = self.values[self.valid_mask]
        if isinstance(self.pixel_area_cm2, np.ndarray):
            masked_area = self.pixel_area_cm2[self.valid_mask]
        else:
            masked_area = float(self.pixel_area_cm2)
        return float(np.sum(masked_values * masked_area))

    def pixel_count_valid(self) -> int:
        """Number of pixels where ``valid_mask`` is ``True``."""
        return int(np.count_nonzero(self.valid_mask))


# ---------------------------------------------------------------------------
# MapGeometry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MapGeometry:
    """Spatial geometry of a quantitative map grid.

    .. note::
        Pixel area is NOT constant for fan-beam DXA scanners because beam
        divergence makes peripheral pixels subtend a larger solid angle than
        central pixels. Use per-pixel ``pixel_area_cm2`` arrays (from the
        calibration step) whenever they are available. The scalar helper
        :meth:`pixel_area_cm2_scalar` is valid only for flat-detector (pencil-
        beam) geometry.
    """

    rows: int
    columns: int
    pixel_spacing_mm: tuple[float, float]
    origin_mm: tuple[float, float]
    is_fan_beam_corrected: bool

    def __post_init__(self) -> None:
        if self.rows <= 0:
            raise ValueError(f"rows must be positive, got {self.rows}.")
        if self.columns <= 0:
            raise ValueError(f"columns must be positive, got {self.columns}.")
        row_sp, col_sp = self.pixel_spacing_mm
        if row_sp <= 0 or col_sp <= 0:
            raise ValueError(
                f"pixel_spacing_mm values must be positive, got {self.pixel_spacing_mm}."
            )

    def pixel_area_cm2_scalar(self) -> float:
        """Approximate scalar pixel area in cm².

        Valid only for flat-detector (pencil-beam) geometry. For fan-beam
        corrected maps use the per-pixel area array from the calibration step.
        """
        row_sp_cm = self.pixel_spacing_mm[0] / 10.0
        col_sp_cm = self.pixel_spacing_mm[1] / 10.0
        return row_sp_cm * col_sp_cm


# ---------------------------------------------------------------------------
# QuantitativeMapBundle
# ---------------------------------------------------------------------------

@dataclass
class QuantitativeMapBundle:
    """Collection of co-registered quantitative maps for one scan.

    At least one of *fat*, *lean*, *bmc*, or *total* must be non-``None``.
    """

    fat: QuantitativeMap | None
    lean: QuantitativeMap | None
    bmc: QuantitativeMap | None
    total: QuantitativeMap | None
    body_mask: np.ndarray | None
    bone_mask: np.ndarray | None
    valid_mask: np.ndarray
    geometry: MapGeometry
    metadata: dict[str, Any]
    provenance: Provenance

    def __post_init__(self) -> None:
        present = [m for m in (self.fat, self.lean, self.bmc, self.total) if m is not None]
        if not present:
            raise ValueError(
                "QuantitativeMapBundle requires at least one non-None map "
                "(fat, lean, bmc, or total)."
            )

        expected_shape = (self.geometry.rows, self.geometry.columns)

        for qmap in present:
            if qmap.values.shape != expected_shape:
                raise ValueError(
                    f"Map '{qmap.quantity}' has shape {qmap.values.shape} but "
                    f"geometry specifies {expected_shape}."
                )

        if self.valid_mask.shape != expected_shape:
            raise ValueError(
                f"valid_mask shape {self.valid_mask.shape} does not match "
                f"geometry shape {expected_shape}."
            )
        if self.valid_mask.dtype != bool:
            object.__setattr__(self, "valid_mask", self.valid_mask.astype(bool))

        if self.body_mask is not None and self.body_mask.shape != expected_shape:
            raise ValueError(
                f"body_mask shape {self.body_mask.shape} does not match "
                f"geometry shape {expected_shape}."
            )
        if self.bone_mask is not None and self.bone_mask.shape != expected_shape:
            raise ValueError(
                f"bone_mask shape {self.bone_mask.shape} does not match "
                f"geometry shape {expected_shape}."
            )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_bundle(self) -> list[str]:
        """Check bundle consistency and return a list of warning strings.

        An empty list means the bundle passed all checks. Callers should
        treat non-empty lists as actionable quality warnings, not hard errors.
        """
        warnings: list[str] = []
        expected_shape = (self.geometry.rows, self.geometry.columns)

        # 1. Shape consistency (defensive — __post_init__ already checked, but
        #    arrays can be mutated after construction).
        for label, qmap in (
            ("fat", self.fat),
            ("lean", self.lean),
            ("bmc", self.bmc),
            ("total", self.total),
        ):
            if qmap is None:
                continue
            if qmap.values.shape != expected_shape:
                warnings.append(
                    f"Shape mismatch: map '{label}' has shape {qmap.values.shape}, "
                    f"expected {expected_shape}."
                )

        # 2. Non-negative values in the valid region.
        for label, qmap in (
            ("fat", self.fat),
            ("lean", self.lean),
            ("bmc", self.bmc),
            ("total", self.total),
        ):
            if qmap is None:
                continue
            combined_mask = qmap.valid_mask & self.valid_mask
            if combined_mask.any():
                min_val = float(qmap.values[combined_mask].min())
                if min_val < 0.0:
                    warnings.append(
                        f"Map '{label}' has negative values in the valid region "
                        f"(min={min_val:.6g}). Physical density cannot be negative."
                    )

        # 3. total ≈ fat + lean + bmc within 1e-6 relative tolerance.
        if self.total is not None and any(
            m is not None for m in (self.fat, self.lean, self.bmc)
        ):
            reconstructed = np.zeros(expected_shape, dtype=np.float64)
            if self.fat is not None:
                reconstructed += self.fat.values
            if self.lean is not None:
                reconstructed += self.lean.values
            if self.bmc is not None:
                reconstructed += self.bmc.values

            combined_mask = self.total.valid_mask & self.valid_mask
            if combined_mask.any():
                total_vals = self.total.values[combined_mask]
                recon_vals = reconstructed[combined_mask]
                denom = np.abs(total_vals)
                denom = np.where(denom < 1e-12, 1e-12, denom)
                rel_err = np.abs(total_vals - recon_vals) / denom
                max_rel_err = float(rel_err.max())
                if max_rel_err > 1e-6:
                    warnings.append(
                        f"Total map deviates from fat+lean+bmc by up to "
                        f"{max_rel_err:.2e} relative error (tolerance 1e-6). "
                        "Check calibration or composition assumption."
                    )

        # 4. Pixel area available if mass calculation might be needed.
        for label, qmap in (
            ("fat", self.fat),
            ("lean", self.lean),
            ("bmc", self.bmc),
            ("total", self.total),
        ):
            if qmap is None:
                continue
            if qmap.pixel_area_cm2 is None:
                warnings.append(
                    f"Map '{label}' has no pixel_area_cm2; "
                    "integrated mass cannot be computed."
                )

        return warnings
