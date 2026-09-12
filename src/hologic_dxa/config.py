"""Runtime configuration loaded from environment variables and optional config files."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, Field, field_validator


class PipelineConfig(BaseModel):
    """All configurable paths and thresholds — no patient data ever stored here."""

    dicom_input_dir: Path | None = Field(
        default=None,
        description="Root directory for recursive DICOM search.",
    )
    dicom_sr_dir: Path | None = Field(
        default=None,
        description="Optional directory containing Hologic Structured Report DICOMs.",
    )
    apex_csv_path: Path | None = Field(
        default=None,
        description="Optional CSV export from Hologic APEX software.",
    )
    hologic_sdk_path: Path | None = Field(
        default=None,
        description="Optional path to Hologic SDK library.",
    )
    hologic_map_export_dir: Path | None = Field(
        default=None,
        description="Optional directory with pre-calibrated quantitative map arrays.",
    )
    bodyloop_point_cloud: Path | None = Field(
        default=None,
        description="Optional BodyLoop mesh or point cloud file.",
    )
    output_dir: Path = Field(
        default=Path("output"),
        description="Base output directory for all pipeline results.",
    )

    # Validation thresholds (provisional — adjust empirically with phantom data)
    mass_tolerance_relative: Annotated[float, Field(gt=0, lt=1)] = Field(
        default=0.05,
        description=(
            "Provisional maximum acceptable relative error between pixel-integrated "
            "mass and APEX regional value. NOT a validated scientific criterion."
        ),
    )
    mass_tolerance_absolute_g: Annotated[float, Field(gt=0)] = Field(
        default=50.0,
        description="Provisional absolute mass tolerance in grams.",
    )

    @field_validator("dicom_input_dir", "dicom_sr_dir", "apex_csv_path",
                     "hologic_sdk_path", "hologic_map_export_dir",
                     "bodyloop_point_cloud", mode="before")
    @classmethod
    def _none_if_empty(cls, v: object) -> object:
        if v == "" or v == "None":
            return None
        return v

    @classmethod
    def from_env(cls) -> "PipelineConfig":
        """Build config from environment variables (DICOM_INPUT_DIR, etc.)."""
        return cls(
            dicom_input_dir=os.environ.get("DICOM_INPUT_DIR"),
            dicom_sr_dir=os.environ.get("DICOM_SR_DIR"),
            apex_csv_path=os.environ.get("APEX_CSV_PATH"),
            hologic_sdk_path=os.environ.get("HOLOGIC_SDK_PATH"),
            hologic_map_export_dir=os.environ.get("HOLOGIC_MAP_EXPORT_DIR"),
            bodyloop_point_cloud=os.environ.get("BODYLOOP_POINT_CLOUD"),
            output_dir=os.environ.get("OUTPUT_DIR", "output"),
        )
