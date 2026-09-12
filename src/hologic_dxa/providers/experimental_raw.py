"""ExperimentalPRProvider — permanently blocked until undocumented parameters are known.

Two energy measurements cannot unambiguously decompose three tissue components
(fat, lean, BMC) without the additional constraints listed in _REQUIRED_PARAMETERS.
See docs/limitations.md#phase-6 for the mathematical argument.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import structlog

from hologic_dxa.providers.base import QuantitativeDataUnavailableError

if TYPE_CHECKING:
    from hologic_dxa.maps.bundle import QuantitativeMapBundle

log = structlog.get_logger(__name__)


class ExperimentalPRProvider:
    name = "experimental_pr"
    experimental = True

    _REQUIRED_PARAMETERS: list[str] = [
        "native_dimensions",
        "numeric_dtype",
        "endianness",
        "low_energy_image",
        "high_energy_image",
        "air_reference_low_energy",
        "air_reference_high_energy",
        "dark_current_low",
        "dark_current_high",
        "detector_gain_map",
        "bad_pixel_mask",
        "table_correction_factors",
        "scatter_correction_method",
        "fan_beam_geometry",
        "fat_lean_calibration_function",
        "bmc_calibration_function",
        "final_units",
    ]

    def can_read(self, source: Path) -> bool:
        return False

    def load(self, source: Path) -> "QuantitativeMapBundle":
        raise QuantitativeDataUnavailableError(
            "ExperimentalPRProvider cannot produce quantitative maps. "
            "The following parameters are required but not yet documented "
            "for any Hologic device: "
            + ", ".join(self._REQUIRED_PARAMETERS)
            + ". See docs/limitations.md and Phase 6 architecture notes."
        )

    def describe_capabilities(self) -> dict[str, object]:
        return {
            "status": "experimental_blocked",
            "experimental": True,
            "missing_parameters": self._REQUIRED_PARAMETERS,
            "note": (
                "Two energy measurements cannot unambiguously decompose three "
                "tissue components (fat, lean, BMC) without additional constraints. "
                "See docs/limitations.md#phase-6."
            ),
        }
