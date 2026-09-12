"""Unit helpers for areal density values.

Uses pint for dimensional analysis so that conversions are explicit and
auditable.  Only units that appear in Hologic APEX exports are recognised;
unknown unit strings are rejected at the boundary rather than silently
propagated.
"""

from __future__ import annotations

import numpy as np
import pint

from hologic_dxa.logging import get_logger

log = get_logger(__name__)

UREG = pint.UnitRegistry()

KNOWN_AREAL_DENSITY_UNITS: frozenset[str] = frozenset(
    {"g/cm2", "g/cm²", "mg/cm2", "mg/cm²"}
)

# Mapping from user-facing unit strings to pint-compatible expressions.
_TO_PINT: dict[str, str] = {
    "g/cm2": "g/cm**2",
    "g/cm²": "g/cm**2",
    "mg/cm2": "mg/cm**2",
    "mg/cm²": "mg/cm**2",
}

_TARGET_PINT = "g/cm**2"


def validate_areal_density_unit(unit: str) -> None:
    """Raise ValueError if unit is not a recognised areal density unit.

    Parameters
    ----------
    unit:
        Unit string to validate, e.g. ``"g/cm2"`` or ``"mg/cm²"``.

    Raises
    ------
    ValueError
        If the unit is not in :data:`KNOWN_AREAL_DENSITY_UNITS`.
    """
    if unit not in KNOWN_AREAL_DENSITY_UNITS:
        raise ValueError(
            f"Unrecognised areal density unit '{unit}'. "
            f"Accepted values: {sorted(KNOWN_AREAL_DENSITY_UNITS)}. "
            "Hologic APEX typically exports in 'g/cm2' or 'mg/cm2'."
        )


def convert_to_g_cm2(
    value: float | np.ndarray, from_unit: str
) -> float | np.ndarray:
    """Convert areal density to g/cm².

    Parameters
    ----------
    value:
        Scalar or array of areal density values.
    from_unit:
        Source unit string (must be in :data:`KNOWN_AREAL_DENSITY_UNITS`).

    Returns
    -------
    float | np.ndarray
        Value(s) expressed in g/cm², same numeric type as input.

    Raises
    ------
    ValueError
        If *from_unit* is not recognised.
    """
    validate_areal_density_unit(from_unit)
    pint_str = _TO_PINT[from_unit]
    qty: pint.Quantity = UREG.Quantity(value, pint_str)
    result: pint.Quantity = qty.to(_TARGET_PINT)
    magnitude: float | np.ndarray = result.magnitude
    log.debug("unit_conversion", from_unit=from_unit, to_unit="g/cm2")
    return magnitude
