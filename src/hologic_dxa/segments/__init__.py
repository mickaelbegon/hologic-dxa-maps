"""Whole-body segment slicing and density profiles from Hologic DXA scans (experimental)."""

from hologic_dxa.segments.apex_scaling import APEX_REGIONS, RegionScaling, rescale_to_apex
from hologic_dxa.segments.bundle_export import bsp_segment_name, to_calibration_bundle
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
    APEX_SCALED_STATUS,
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
    "APEX_REGIONS",
    "APEX_SCALED_STATUS",
    "CALIBRATION_STATUS",
    "MARKER_NAMES",
    "SEGMENT_DEFINITIONS",
    "XIPHOID_FRACTION",
    "RegionScaling",
    "SegmentDefinition",
    "SegmentProfile",
    "SliceGeometry",
    "SliceProfile",
    "TissueModel",
    "WholeBodyMu",
    "bsp_segment_name",
    "build_segments",
    "compute_segment_profiles",
    "lateral_constraint",
    "load_cor_json",
    "load_whole_body_mu",
    "parse_r_file_images",
    "profiles_to_dataframe",
    "rescale_to_apex",
    "slice_segment",
    "suggest_markers",
    "to_calibration_bundle",
]
