"""Per-pixel quality control checks for quantitative DXA map bundles.

All functions return a list of warning strings.  An empty list means the
check passed.  Callers decide what to do with warnings — log them, surface
them in a validation report, or raise.
"""

from __future__ import annotations

import numpy as np

from hologic_dxa.logging import get_logger

log = get_logger(__name__)


def check_total_consistency(
    fat: np.ndarray,
    lean: np.ndarray,
    bmc: np.ndarray,
    total: np.ndarray,
    valid_mask: np.ndarray,
    rtol: float = 1e-5,
) -> list[str]:
    """Return list of QC warnings.  Empty list means all checks passed.

    Checks performed:
    1. ``total ≈ fat + lean + bmc`` in every valid pixel (relative tolerance).
    2. No negative values in fat, lean, bmc, or total within the valid region.

    Parameters
    ----------
    fat, lean, bmc, total:
        2-D areal density arrays (g/cm²), same shape.
    valid_mask:
        Boolean array; True pixels are included in the checks.
    rtol:
        Maximum acceptable relative deviation for the consistency check.

    Returns
    -------
    list[str]
        Warning messages.  Empty when all checks pass.
    """
    warnings: list[str] = []

    computed_total = fat + lean + bmc
    abs_diff = np.abs(total - computed_total)
    denominator = np.maximum(np.abs(total), 1e-12)
    rel_diff = abs_diff / denominator

    inconsistent = valid_mask & (rel_diff > rtol)
    if np.any(inconsistent):
        n_bad = int(np.sum(inconsistent))
        max_rel = float(np.max(rel_diff[inconsistent]))
        msg = (
            f"total ≠ fat + lean + bmc in {n_bad} valid pixel(s); "
            f"maximum relative deviation {max_rel:.3e} (rtol={rtol:.1e}). "
            "This may indicate a calibration inconsistency or a rounding artefact."
        )
        warnings.append(msg)
        log.warning("qc_total_consistency_failed", n_pixels=n_bad, max_rel_dev=max_rel)

    for arr, name in ((fat, "fat"), (lean, "lean"), (bmc, "bmc"), (total, "total")):
        warnings.extend(check_non_negative(arr, name, valid_mask))

    return warnings


def check_non_negative(
    arr: np.ndarray,
    name: str,
    valid_mask: np.ndarray,
    atol: float = 1e-6,
) -> list[str]:
    """Return warnings for pixels with value < -atol in the valid region.

    Parameters
    ----------
    arr:
        2-D array to check.
    name:
        Human-readable name used in the warning message.
    valid_mask:
        Boolean inclusion mask.
    atol:
        Absolute tolerance; values in ``(-atol, 0]`` are accepted without warning.

    Returns
    -------
    list[str]
        Warning messages.  Empty when all valid pixels are ≥ -atol.
    """
    negative = valid_mask & (arr < -atol)
    if not np.any(negative):
        return []

    n_neg = int(np.sum(negative))
    min_val = float(np.min(arr[negative]))
    msg = (
        f"'{name}' contains {n_neg} valid pixel(s) with value < -{atol:.1e}; "
        f"minimum value: {min_val:.6g} g/cm². "
        "Negative areal density values are physically meaningless and indicate "
        "a calibration or reconstruction error."
    )
    log.warning("qc_negative_values", array=name, n_pixels=n_neg, min_value=min_val)
    return [msg]
