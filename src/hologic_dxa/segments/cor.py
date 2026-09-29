"""Centres of rotation (CoR): marker names, segment definitions, JSON loading."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

MARKER_NAMES: tuple[str, ...] = (
    "head_top",
    "atlanto_occipital",
    "cervicothoracic",
    "xiphoid",
    "umbilicus",
    "lumbosacral",
    "crotch",
    "glenohumeral_L",
    "elbow_L",
    "wrist_L",
    "hand_end_L",
    "hip_L",
    "knee_L",
    "ankle_L",
    "foot_end_L",
    "glenohumeral_R",
    "elbow_R",
    "wrist_R",
    "hand_end_R",
    "hip_R",
    "knee_R",
    "ankle_R",
    "foot_end_R",
)

# xiphoid, umbilicus and crotch are trunk landmarks reserved for a future thorax /
# abdomen split (Hatze, Yeadon); no segment uses them yet.
# Markers derived from others (never read from JSON).
VIRTUAL_MARKERS: dict[str, tuple[str, str]] = {"hip_mid": ("hip_L", "hip_R")}

Cor = dict[str, tuple[float, float]]


@dataclass(frozen=True)
class SegmentDefinition:
    name: str
    proximal: str
    distal: str
    tissue: str
    uniform_density: bool = False


SEGMENT_DEFINITIONS: tuple[SegmentDefinition, ...] = (
    SegmentDefinition("head", "head_top", "atlanto_occipital", "head"),
    SegmentDefinition("neck", "atlanto_occipital", "cervicothoracic", "neck"),
    SegmentDefinition("trunk", "cervicothoracic", "lumbosacral", "trunk"),
    SegmentDefinition("pelvis", "lumbosacral", "hip_mid", "pelvis"),
    SegmentDefinition("upper_arm_L", "glenohumeral_L", "elbow_L", "upper_arm"),
    SegmentDefinition("upper_arm_R", "glenohumeral_R", "elbow_R", "upper_arm"),
    SegmentDefinition("forearm_L", "elbow_L", "wrist_L", "forearm"),
    SegmentDefinition("forearm_R", "elbow_R", "wrist_R", "forearm"),
    SegmentDefinition("hand_L", "wrist_L", "hand_end_L", "hand", uniform_density=True),
    SegmentDefinition("hand_R", "wrist_R", "hand_end_R", "hand", uniform_density=True),
    SegmentDefinition("thigh_L", "hip_L", "knee_L", "thigh"),
    SegmentDefinition("thigh_R", "hip_R", "knee_R", "thigh"),
    SegmentDefinition("shank_L", "knee_L", "ankle_L", "shank"),
    SegmentDefinition("shank_R", "knee_R", "ankle_R", "shank"),
    SegmentDefinition("foot_L", "ankle_L", "foot_end_L", "foot", uniform_density=True),
    SegmentDefinition("foot_R", "ankle_R", "foot_end_R", "foot", uniform_density=True),
)


def load_cor_json(path: Path, grid_shape: tuple[int, int] = (106, 150)) -> Cor:
    """Load ``{marker: [row, col]}`` in whole-body grid pixels (row 0 = patient right).

    Raises
    ------
    ValueError
        On unknown marker names, malformed values, or points outside the grid.
    """
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: expected a JSON object {{marker: [row, col]}}.")
    cor: Cor = {}
    for name, value in raw.items():
        if name not in MARKER_NAMES:
            raise ValueError(f"{path}: unknown marker '{name}'. Known: {', '.join(MARKER_NAMES)}.")
        if not (isinstance(value, list) and len(value) == 2):
            raise ValueError(f"{path}: marker '{name}' must be [row, col].")
        row, col = float(value[0]), float(value[1])
        if not (math.isfinite(row) and math.isfinite(col)):
            raise ValueError(f"{path}: marker '{name}' has a non-finite coordinate.")
        if not (0 <= row < grid_shape[0] and 0 <= col < grid_shape[1]):
            raise ValueError(
                f"{path}: marker '{name}' = ({row}, {col}) is outside the "
                f"{grid_shape[0]}x{grid_shape[1]} grid."
            )
        cor[name] = (row, col)
    return cor


def with_virtual_markers(cor: Cor) -> Cor:
    """Return a copy of ``cor`` extended with derived markers (e.g. ``hip_mid``)."""
    out = dict(cor)
    for name, (a, b) in VIRTUAL_MARKERS.items():
        if a in cor and b in cor:
            out[name] = (0.5 * (cor[a][0] + cor[b][0]), 0.5 * (cor[a][1] + cor[b][1]))
    return out


def build_segments(
    cor: Cor,
) -> list[tuple[SegmentDefinition, tuple[float, float], tuple[float, float]]]:
    """Return (definition, proximal_px, distal_px) for every segment whose two CoR exist."""
    full = with_virtual_markers(cor)
    return [
        (d, full[d.proximal], full[d.distal])
        for d in SEGMENT_DEFINITIONS
        if d.proximal in full and d.distal in full
    ]
