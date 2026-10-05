"""Export segment density profiles as a BodyLoop ``DXACalibrationBundle`` JSON.

The bridge is a plain JSON file, so the consumer needs no dependency on this
package. Keys follow ``bodyloop_anthropometrics.anthropometry.dxa_calibration``:
``regional_fat_fraction`` (Tier 1), ``segment_lengths_m`` and ``slice_fat_fraction``
(Tier 3, proximal to distal).

Nothing is imputed: a segment with any invalid slice is left out of
``slice_fat_fraction`` and listed, with its reason, in ``omitted_segments``.
``regional_bmc_g`` is not produced (no calibrated bone map here), so Tier 2
cannot be run from this bundle.
"""

from __future__ import annotations

from typing import Any

from hologic_dxa.segments.profiles import APEX_SCALED_STATUS, CALIBRATION_STATUS, SegmentProfile

LIMB_TYPES = ("upper_arm", "forearm", "thigh", "shank")
SEXES = ("male", "female", "other")


def bsp_segment_name(name: str) -> str:
    """``thigh_L`` -> ``thigh_left``; names without a side are unchanged."""
    if name.endswith("_L"):
        return name[:-2] + "_left"
    if name.endswith("_R"):
        return name[:-2] + "_right"
    return name


def to_calibration_bundle(
    profiles: list[SegmentProfile],
    *,
    sex: str,
    provenance: str,
    n_slices: int,
) -> dict[str, Any]:
    """Build the bundle dictionary from segment profiles.

    ``regional_fat_fraction`` pools left and right of each limb type, weighted by
    the ``mu_H`` integral of each slice (an attenuation proxy for mass); it is a
    soft-tissue fat fraction from the empirical R relation, not a Hologic value.

    Raises
    ------
    ValueError
        On an unknown ``sex`` or when no limb segment is available.
    """
    if sex not in SEXES:
        raise ValueError(f"sex must be one of {SEXES}, got {sex!r}")

    slice_fat: dict[str, list[float]] = {}
    lengths: dict[str, float] = {}
    omitted: dict[str, str] = {}
    pooled: dict[str, list[float]] = {t: [0.0, 0.0] for t in LIMB_TYPES}  # sum(f*w), sum(w)

    for seg in profiles:
        bsp_name = bsp_segment_name(seg.name)
        lengths[bsp_name] = seg.length_m
        valid = [s for s in seg.slices if s is not None]
        uniform = bool(valid) and all(s.density_source == "uniform_segment_mean" for s in valid)
        if not valid:
            omitted[bsp_name] = "no valid slice"
        elif uniform:
            slice_fat[bsp_name] = [valid[0].fat_fraction] * n_slices
        elif len(valid) < len(seg.slices):
            n_bad = len(seg.slices) - len(valid)
            omitted[bsp_name] = f"{n_bad} of {len(seg.slices)} slices invalid"
        else:
            slice_fat[bsp_name] = [s.fat_fraction for s in valid]
        if seg.tissue in pooled:
            for s in valid:
                pooled[seg.tissue][0] += s.fat_fraction * s.mu_h_integral
                pooled[seg.tissue][1] += s.mu_h_integral

    regional = {t: num / den for t, (num, den) in pooled.items() if den > 0}
    if not regional:
        raise ValueError("No limb segment (upper_arm, forearm, thigh, shank) has valid slices.")

    # The status describes the limb segments BodyLoop consumes; head, neck, trunk and
    # pelvis are never rescaled to APEX.
    exported = [
        s
        for seg in profiles
        if seg.tissue in (*LIMB_TYPES, "hand", "foot") and bsp_segment_name(seg.name) in slice_fat
        for s in seg.slices
        if s is not None
    ]
    all_scaled = bool(exported) and all(s.scaled_to_apex for s in exported)

    return {
        "sex": sex,
        "provenance": provenance,
        "calibration_status": APEX_SCALED_STATUS if all_scaled else CALIBRATION_STATUS,
        "slice_orientation": "proximal_to_distal",
        "regional_fat_fraction": regional,
        "segment_lengths_m": lengths,
        "slice_fat_fraction": slice_fat,
        "omitted_segments": omitted,
    }
