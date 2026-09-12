"""QuantitativeMapBundle — canonical container for calibrated spatial tissue maps."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class QuantitativeMapBundle:
    """Holds calibrated, geometry-corrected spatial maps for all tissue components.

    All arrays must share the same shape.  Values are in g/cm².
    `total` is a derived property and is never stored separately.
    """

    fat: np.ndarray
    lean: np.ndarray
    bmc: np.ndarray
    valid_mask: np.ndarray
    pixel_area_cm2: float
    source_sop_instance_uid: str
    already_geometry_corrected: bool = True
    provider_name: str = ""
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        shapes = {arr.shape for arr in (self.fat, self.lean, self.bmc, self.valid_mask)}
        if len(shapes) != 1:
            raise ValueError(
                f"All component arrays must share one shape; got {shapes}"
            )

    @property
    def total(self) -> np.ndarray:
        """Derived total mass density (fat + lean + bmc), never stored separately."""
        return self.fat + self.lean + self.bmc

    @property
    def shape(self) -> tuple[int, ...]:
        return self.fat.shape  # type: ignore[return-value]
