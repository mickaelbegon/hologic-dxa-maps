"""Phantom calibration interface for bi-energy DXA reconstruction.

A calibration phantom with known fat-equivalent and bone-equivalent
material masses provides the calibration function:

    (a_L, a_H) = F_θ(σ_fat, σ_lean)   for soft tissue pixels
    (a_L, a_H) = F_θ(σ_fat, σ_lean, σ_BMC)  for bone pixels

where a_L and a_H are the log-attenuation values at low and high energies.

The calibration is device-specific, APEX-version-specific, and potentially
patient-size-dependent. It MUST be derived from acquisitions on calibration
phantoms with certified composition, not from a single human scan.

This module is BLOCKED until phantom data with matching archive DICOMs
(P/R files) and certified reference values are available.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class PhantomCalibrationPoint:
    """One point in the calibration space."""

    sigma_fat_g_cm2: float
    sigma_lean_g_cm2: float
    sigma_bmc_g_cm2: float      # 0 for soft-tissue-only calibration
    attenuation_low: float      # a_L = -ln(I_L / I_0L)
    attenuation_high: float     # a_H = -ln(I_H / I_0H)
    phantom_id: str
    notes: str = ""


class CalibrationFunction:
    """Fitted bi-energy calibration function.

    F_θ: (σ_fat, σ_lean) → (a_L, a_H) for soft tissue
    F_θ: (σ_fat, σ_lean, σ_BMC) → (a_L, a_H) for bone (underdetermined)

    Fitting requires at minimum:
    - 4 soft-tissue calibration points (fat equivalent, lean equivalent,
      and two mixtures) for the 2D → 2D mapping
    - Additional bone calibration points for the 3D case

    Status: NOT IMPLEMENTED — requires phantom acquisition data.
    """

    def __init__(self) -> None:
        self._fitted = False

    def fit(self, calibration_points: list[PhantomCalibrationPoint]) -> None:
        raise NotImplementedError(
            "Calibration fitting is not implemented. "
            "Required: phantom acquisition DICOMs with P/R files and "
            "certified reference composition values. "
            "Minimum: 4 soft-tissue calibration points for 2×2 inversion, "
            "plus additional bone points for the underdetermined 3-component case."
        )

    def predict(
        self,
        sigma_fat: np.ndarray,
        sigma_lean: np.ndarray,
        sigma_bmc: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        if not self._fitted:
            raise RuntimeError(
                "CalibrationFunction has not been fitted. Call fit() first."
            )
        raise NotImplementedError("fit() was not completed.")

    @classmethod
    def load(cls, path: Path) -> "CalibrationFunction":
        raise NotImplementedError(
            "Calibration serialization is not yet implemented."
        )
