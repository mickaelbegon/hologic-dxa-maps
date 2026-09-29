"""Whole-body segment slicing and density profiles from Hologic DXA scans (experimental)."""

from hologic_dxa.segments.cor import (
    MARKER_NAMES,
    SEGMENT_DEFINITIONS,
    XIPHOID_FRACTION,
    SegmentDefinition,
    build_segments,
    load_cor_json,
    suggest_markers,
)
from hologic_dxa.segments.profiles import (
    CALIBRATION_STATUS,
    SegmentProfile,
    SliceGeometry,
    SliceProfile,
    TissueModel,
    compute_segment_profiles,
    lateral_constraint,
    profiles_to_dataframe,
    slice_segment,
)
from hologic_dxa.segments.whole_body import (
    WholeBodyMu,
    load_whole_body_mu,
    parse_r_file_images,
)

__all__ = [
    "CALIBRATION_STATUS",
    "MARKER_NAMES",
    "SEGMENT_DEFINITIONS",
    "XIPHOID_FRACTION",
    "SegmentDefinition",
    "SegmentProfile",
    "SliceGeometry",
    "SliceProfile",
    "TissueModel",
    "WholeBodyMu",
    "build_segments",
    "compute_segment_profiles",
    "lateral_constraint",
    "load_cor_json",
    "load_whole_body_mu",
    "parse_r_file_images",
    "profiles_to_dataframe",
    "slice_segment",
    "suggest_markers",
]
