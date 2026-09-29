"""Perpendicular slicing of body segments and 2-component density profiles.

EXPERIMENTAL / NOT CALIBRATED: fat fraction comes from an empirical linear
relation on R = mu_L / mu_H, and the bone fractions and tissue densities are
literature-style assumptions (see ``TissueModel``). Outputs are relative
density profiles, not validated absolute densities.

Physical space: x = lateral (row * row_size_m, row 0 = patient right),
y = head-to-foot (col * col_size_m, col 0 = top of head).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates

from hologic_dxa.segments.cor import Cor, SegmentDefinition, build_segments, with_virtual_markers

CALIBRATION_STATUS = "empirical_uncalibrated"

AllowedFn = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True)
class SliceGeometry:
    row_size_m: float = 0.0057
    col_size_m: float = 0.0131
    body_threshold: float = 50.0
    n_samples: int = 300
    half_window_m: float = 0.35
    # Max half-width [m] of a slice around the segment axis, per tissue type.
    half_width_cap_m: Mapping[str, float] = field(
        default_factory=lambda: {
            "head": 0.10,
            "neck": 0.07,
            "trunk": 0.30,
            "pelvis": 0.22,
            "upper_arm": 0.07,
            "forearm": 0.06,
            "hand": 0.05,
            "thigh": 0.12,
            "shank": 0.09,
            "foot": 0.06,
        }
    )
    # Tissues whose medial extent from the axis is limited to their lateral extent
    # (prevents an arm held against the trunk from absorbing trunk pixels).
    symmetric_edge_tissues: frozenset[str] = frozenset({"upper_arm"})
    # Pelvis slices stay within the hip-CoR half-span plus this margin [m].
    pelvis_margin_m: float = 0.09


@dataclass(frozen=True)
class TissueModel:
    cal_slope: float = -2.49
    cal_intercept: float = 3.41
    rho_fat: float = 916.0
    rho_muscle: float = 1056.0
    rho_bone: float = 1560.0
    bone_fraction: Mapping[str, float] = field(
        default_factory=lambda: {
            "head": 0.095,
            "neck": 0.06,
            "trunk": 0.065,
            "pelvis": 0.08,
            "upper_arm": 0.045,
            "forearm": 0.055,
            "hand": 0.10,
            "thigh": 0.060,
            "shank": 0.055,
            "foot": 0.10,
        }
    )
    default_bone_fraction: float = 0.06
    # tissue -> fraction of proximal slices whose density is taken from the distal side.
    borrow_distal_density: Mapping[str, float] = field(
        default_factory=lambda: {"upper_arm": 0.5}
    )

    def density(self, r_mean: float, tissue: str) -> tuple[float, float]:
        """Return (fat_fraction, density_kg_m3) for a mean R = mu_L / mu_H."""
        f_bone = self.bone_fraction.get(tissue, self.default_bone_fraction)
        f_fat = float(np.clip(self.cal_slope * r_mean + self.cal_intercept, 0.0, 1.0))
        f_soft = 1.0 - f_bone
        rho = (
            self.rho_bone * f_bone
            + self.rho_fat * f_fat * f_soft
            + self.rho_muscle * (1.0 - f_fat) * f_soft
        )
        return f_fat, float(rho)


@dataclass(frozen=True)
class SliceProfile:
    slice_idx: int
    height_along_seg_m: float
    center_row: float
    center_col: float
    edge_rows: tuple[float, float]
    edge_cols: tuple[float, float]
    n_body: int
    body_width_m: float
    r_mean: float
    fat_fraction: float
    mu_h_integral: float
    density_kg_m3: float
    density_source: str = "measured"


@dataclass(frozen=True)
class SegmentProfile:
    name: str
    tissue: str
    proximal: tuple[float, float]
    distal: tuple[float, float]
    length_m: float
    slices: list[SliceProfile | None]


def _median_axis_x(cor: Cor, geometry: SliceGeometry) -> Callable[[np.ndarray], np.ndarray] | None:
    """x_mid(y): median axis through C7 and L5/S1, or the hip midpoint as fallback."""
    if "cervicothoracic" in cor and "lumbosacral" in cor:
        (r0, c0), (r1, c1) = cor["cervicothoracic"], cor["lumbosacral"]
        x0, y0 = r0 * geometry.row_size_m, c0 * geometry.col_size_m
        x1, y1 = r1 * geometry.row_size_m, c1 * geometry.col_size_m
        if abs(y1 - y0) > 1e-6:
            slope = (x1 - x0) / (y1 - y0)
            return lambda y: x0 + (y - y0) * slope
        return lambda y: np.full_like(y, 0.5 * (x0 + x1))
    if "hip_L" in cor and "hip_R" in cor:
        x_mid = 0.5 * (cor["hip_L"][0] + cor["hip_R"][0]) * geometry.row_size_m
        return lambda y: np.full_like(y, x_mid)
    return None


def lateral_constraint(
    definition: SegmentDefinition, cor: Cor, geometry: SliceGeometry
) -> AllowedFn | None:
    """Anatomical lateral limit for a segment, as a mask function of physical (x, y).

    - trunk: between the two glenohumeral CoR (no spill onto the arms);
    - pelvis: hip-CoR half-span plus a margin around the median axis (no spill onto
      the arms or hands resting against the hips);
    - thigh: on its own side of the median axis (no spill onto the other thigh).
    """
    cor = with_virtual_markers(cor)
    if definition.tissue == "trunk":
        if "glenohumeral_L" in cor and "glenohumeral_R" in cor:
            xs = sorted(
                cor[k][0] * geometry.row_size_m for k in ("glenohumeral_L", "glenohumeral_R")
            )
            return lambda x, y: (x >= xs[0]) & (x <= xs[1])
        return None
    if definition.tissue == "pelvis":
        x_mid = _median_axis_x(cor, geometry)
        if x_mid is None or "hip_L" not in cor or "hip_R" not in cor:
            return None
        half = (
            abs(cor["hip_L"][0] - cor["hip_R"][0]) * geometry.row_size_m / 2
            + geometry.pelvis_margin_m
        )
        return lambda x, y: np.abs(x - x_mid(y)) <= half
    if definition.tissue == "thigh":
        x_mid = _median_axis_x(cor, geometry)
        if x_mid is None:
            return None
        sign = 1.0 if definition.name.endswith("_L") else -1.0  # left = larger row
        return lambda x, y: sign * (x - x_mid(y)) > 0
    return None


def slice_segment(
    definition: SegmentDefinition,
    proximal_px: tuple[float, float],
    distal_px: tuple[float, float],
    mu_h: np.ndarray,
    mu_l: np.ndarray,
    n_slices: int = 10,
    geometry: SliceGeometry | None = None,
    tissue_model: TissueModel | None = None,
    allowed: AllowedFn | None = None,
) -> list[SliceProfile | None]:
    """Cut one segment into ``n_slices`` slices perpendicular to its CoR-to-CoR axis."""
    geom = geometry or SliceGeometry()
    model = tissue_model or TissueModel()
    n_rows, n_cols = mu_h.shape
    scale = np.array([geom.row_size_m, geom.col_size_m])
    prox_m = np.asarray(proximal_px, dtype=float) * scale
    axis_m = np.asarray(distal_px, dtype=float) * scale - prox_m
    axis_len = float(np.linalg.norm(axis_m))
    if axis_len < 0.02:
        return [None] * n_slices

    axis_u = axis_m / axis_len
    perp_u = np.array([-axis_u[1], axis_u[0]])
    t_m = np.linspace(-geom.half_window_m, geom.half_window_m, geom.n_samples)
    dt_m = abs(float(t_m[1] - t_m[0]))
    cap = geom.half_width_cap_m.get(definition.tissue, 0.12)
    centre_i = len(t_m) // 2

    out_sign = 0.0
    if (
        definition.tissue in geom.symmetric_edge_tissues
        and definition.name.endswith(("_L", "_R"))
        and abs(perp_u[0]) > 1e-6
    ):
        outward_x = 1.0 if definition.name.endswith("_L") else -1.0  # left = larger row
        out_sign = float(np.sign(perp_u[0])) * outward_x

    out: list[SliceProfile | None] = []
    for i in range(n_slices):
        frac = (i + 0.5) / n_slices
        center_m = prox_m + frac * axis_m
        pts = center_m + np.outer(t_m, perp_u)
        rows, cols = pts[:, 0] / geom.row_size_m, pts[:, 1] / geom.col_size_m

        valid = (rows >= 0) & (rows < n_rows) & (cols >= 0) & (cols < n_cols)
        rc = np.clip(rows, 0, n_rows - 1e-3)
        cc = np.clip(cols, 0, n_cols - 1e-3)
        v_h = map_coordinates(mu_h, [rc, cc], order=1, prefilter=False)
        v_l = map_coordinates(mu_l, [rc, cc], order=1, prefilter=False)
        v_h[~valid] = 0.0
        v_l[~valid] = 0.0

        body = (v_h > geom.body_threshold) & valid & (np.abs(t_m) <= cap)
        if allowed is not None:
            body &= allowed(pts[:, 0], pts[:, 1])
        idx = np.where(body)[0]
        if idx.size < 5:
            out.append(None)
            continue
        runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        best = min(
            runs,
            key=lambda r: 0
            if r[0] <= centre_i <= r[-1]
            else min(abs(int(r[0]) - centre_i), abs(int(r[-1]) - centre_i)),
        )
        if out_sign != 0.0:
            t_signed = t_m[best] * out_sign  # > 0 on the outer (lateral) side
            d_out = float(t_signed.max())
            if d_out > 0:
                best = best[t_signed >= -d_out]
        if best.size < 5:
            out.append(None)
            continue

        vh_b, vl_b = v_h[best], v_l[best]
        r_mean = float(vl_b.sum() / max(vh_b.sum(), 1e-6))
        f_fat, rho = model.density(r_mean, definition.tissue)
        out.append(
            SliceProfile(
                slice_idx=i,
                height_along_seg_m=frac * axis_len,
                center_row=float(center_m[0] / geom.row_size_m),
                center_col=float(center_m[1] / geom.col_size_m),
                edge_rows=(float(rows[best[0]]), float(rows[best[-1]])),
                edge_cols=(float(cols[best[0]]), float(cols[best[-1]])),
                n_body=int(best.size),
                body_width_m=float(best.size * dt_m),
                r_mean=r_mean,
                fat_fraction=f_fat,
                mu_h_integral=float(vh_b.sum() * dt_m),
                density_kg_m3=rho,
            )
        )
    return out


def _borrow_distal_density(
    slices: list[SliceProfile | None], fraction: float
) -> list[SliceProfile | None]:
    """Proximal slices take r_mean/fat/density from the first valid distal-side slice."""
    n_borrow = int(fraction * len(slices))
    ref = next((s for s in slices[n_borrow:] if s is not None), None)
    if ref is None or n_borrow == 0:
        return slices
    return [
        replace(
            s,
            r_mean=ref.r_mean,
            fat_fraction=ref.fat_fraction,
            density_kg_m3=ref.density_kg_m3,
            density_source="borrowed_distal",
        )
        if s is not None and i < n_borrow
        else s
        for i, s in enumerate(slices)
    ]


def _uniform_density(
    slices: list[SliceProfile | None], tissue: str, model: TissueModel
) -> list[SliceProfile | None]:
    """Give every slice the segment-pooled R (weighted by mu_H integral)."""
    valid = [s for s in slices if s is not None]
    weight = sum(s.mu_h_integral for s in valid)
    if not valid or weight <= 0:
        return slices
    r_pool = sum(s.r_mean * s.mu_h_integral for s in valid) / weight
    f_fat, rho = model.density(r_pool, tissue)
    return [
        replace(
            s,
            r_mean=r_pool,
            fat_fraction=f_fat,
            density_kg_m3=rho,
            density_source="uniform_segment_mean",
        )
        if s is not None
        else None
        for s in slices
    ]


def compute_segment_profiles(
    mu_h: np.ndarray,
    mu_l: np.ndarray,
    cor: Cor,
    n_slices: int = 10,
    geometry: SliceGeometry | None = None,
    tissue_model: TissueModel | None = None,
) -> list[SegmentProfile]:
    """Density profiles for every segment whose two bounding CoR are present."""
    if mu_h.shape != mu_l.shape or mu_h.ndim != 2:
        raise ValueError(
            f"mu_h and mu_l must be 2-D and same shape; got {mu_h.shape}, {mu_l.shape}."
        )
    geom = geometry or SliceGeometry()
    model = tissue_model or TissueModel()
    scale = np.array([geom.row_size_m, geom.col_size_m])
    profiles: list[SegmentProfile] = []
    for definition, prox, dist in build_segments(cor):
        length = float(np.linalg.norm((np.array(dist) - np.array(prox)) * scale))
        slices = slice_segment(
            definition,
            prox,
            dist,
            mu_h,
            mu_l,
            n_slices=n_slices,
            geometry=geom,
            tissue_model=model,
            allowed=lateral_constraint(definition, cor, geom),
        )
        fraction = model.borrow_distal_density.get(definition.tissue, 0.0)
        if fraction > 0:
            slices = _borrow_distal_density(slices, fraction)
        if definition.uniform_density:
            slices = _uniform_density(slices, definition.tissue, model)
        profiles.append(
            SegmentProfile(definition.name, definition.tissue, prox, dist, length, slices)
        )
    return profiles


def profiles_to_dataframe(profiles: list[SegmentProfile]) -> pd.DataFrame:
    """One row per valid slice; always carries the calibration-status flag."""
    rows = []
    for seg in profiles:
        for s in seg.slices:
            if s is None:
                continue
            rows.append(
                {
                    "segment": seg.name,
                    "tissue": seg.tissue,
                    "slice_idx": s.slice_idx,
                    "height_along_seg_m": s.height_along_seg_m,
                    "body_width_m": s.body_width_m,
                    "r_mean": s.r_mean,
                    "fat_fraction": s.fat_fraction,
                    "mu_h_integral": s.mu_h_integral,
                    "density_kg_m3": s.density_kg_m3,
                    "density_source": s.density_source,
                    "calibration_status": CALIBRATION_STATUS,
                }
            )
    return pd.DataFrame(rows)
