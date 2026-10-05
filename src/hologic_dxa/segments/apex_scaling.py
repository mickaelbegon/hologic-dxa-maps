"""Rescale limb fat fractions to the Hologic APEX regional values.

The slice fat fractions come from an empirical, uncalibrated relation on
R = mu_L / mu_H. APEX reports calibrated fat and lean mass for each arm and leg;
this module keeps the *shape* of the slice profile along the limb and multiplies it
by one factor per region so that the mass-weighted soft-tissue fat fraction equals
``fat_g / (fat_g + lean_g)``.

Assumptions (documented, not validated):
- the APEX region "L Arm" / "L Leg" is the patient's left limb and covers the hand /
  foot, so it maps to upper arm + forearm + hand / thigh + shank + foot;
- the weight of a slice is its ``mu_H`` integral times the slice thickness
  (segment length / number of slices), a proxy for tissue mass;
- head and trunk are not rescaled: their APEX region boundaries do not match the
  segments here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace

from hologic_dxa.segments.profiles import SegmentProfile, SliceProfile, TissueModel

APEX_REGIONS: dict[str, tuple[str, ...]] = {
    "L Arm": ("upper_arm_L", "forearm_L", "hand_L"),
    "R Arm": ("upper_arm_R", "forearm_R", "hand_R"),
    "L Leg": ("thigh_L", "shank_L", "foot_L"),
    "R Leg": ("thigh_R", "shank_R", "foot_R"),
}

_TOLERANCE = 1e-4


@dataclass(frozen=True)
class RegionScaling:
    region: str
    fat_fraction_slices: float
    fat_fraction_apex: float
    factor: float


def _weighted_fraction(items: list[tuple[float, float]], factor: float) -> float:
    total = sum(w for _, w in items)
    return sum(min(1.0, factor * f) * w for f, w in items) / total


def _solve_factor(items: list[tuple[float, float]], target: float) -> float | None:
    """Factor k with weighted mean of min(1, k f) equal to target, or None if unreachable."""
    hi = 100.0
    if _weighted_fraction(items, hi) < target - _TOLERANCE:
        return None
    lo = 0.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _weighted_fraction(items, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rescale_to_apex(
    profiles: list[SegmentProfile],
    apex_fat_lean_g: Mapping[str, tuple[float, float]],
    tissue_model: TissueModel | None = None,
) -> tuple[list[SegmentProfile], list[RegionScaling], list[str]]:
    """Scale each limb region to its APEX ``fat / (fat + lean)``.

    Parameters
    ----------
    profiles:
        Output of ``compute_segment_profiles``.
    apex_fat_lean_g:
        ``{region: (fat_g, lean_g)}`` for "L Arm", "R Arm", "L Leg", "R Leg".
    tissue_model:
        Used to recompute the density from the scaled fat fraction.

    Returns
    -------
    (profiles, scalings, warnings)
        New profiles (inputs are not modified), one ``RegionScaling`` per scaled
        region, and messages for every region left unscaled and why.
    """
    model = tissue_model or TissueModel()
    by_name = {p.name: p for p in profiles}
    out = dict(by_name)
    scalings: list[RegionScaling] = []
    warnings: list[str] = []

    for region, names in APEX_REGIONS.items():
        missing = [n for n in names if n not in by_name]
        if missing:
            warnings.append(f"{region}: not scaled, segments missing ({', '.join(missing)}).")
            continue
        if region not in apex_fat_lean_g:
            warnings.append(f"{region}: not scaled, no APEX fat/lean values.")
            continue
        fat_g, lean_g = apex_fat_lean_g[region]
        if fat_g is None or lean_g is None or fat_g + lean_g <= 0:
            warnings.append(f"{region}: not scaled, invalid APEX fat/lean values.")
            continue
        target = fat_g / (fat_g + lean_g)

        items: list[tuple[float, float]] = []
        for n in names:
            seg = by_name[n]
            thickness = seg.length_m / len(seg.slices)
            items += [
                (s.fat_fraction, s.mu_h_integral * thickness) for s in seg.slices if s is not None
            ]
        if not items or sum(w for _, w in items) <= 0:
            warnings.append(f"{region}: not scaled, no valid slice.")
            continue
        factor = _solve_factor(items, target)
        if factor is None:
            warnings.append(f"{region}: not scaled, APEX fat fraction {target:.3f} unreachable.")
            continue

        for n in names:
            seg = by_name[n]
            new_slices: list[SliceProfile | None] = []
            for s in seg.slices:
                if s is None:
                    new_slices.append(None)
                    continue
                f_new = min(1.0, factor * s.fat_fraction)
                new_slices.append(
                    replace(
                        s,
                        fat_fraction=f_new,
                        density_kg_m3=model.density_from_fat_fraction(f_new, seg.tissue),
                        fat_scale=factor,
                        scaled_to_apex=True,
                    )
                )
            out[n] = replace(seg, slices=new_slices)
        scalings.append(
            RegionScaling(region, _weighted_fraction(items, 1.0), target, factor)
        )

    return [out[p.name] for p in profiles], scalings, warnings
