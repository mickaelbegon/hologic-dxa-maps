"""HologicSDKProvider — stub pending Hologic SDK acquisition.

This module intentionally raises NotImplementedError from load() because the
Hologic SDK, its license, and its API documentation have not been provided.
The describe_capabilities() method records exactly what is needed to unblock it.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import structlog
from pydantic import BaseModel

if TYPE_CHECKING:
    from hologic_dxa.maps.bundle import QuantitativeMapBundle

log = structlog.get_logger(__name__)


class HologicSDKConfig(BaseModel):
    sdk_path: Path
    license_file: Path | None = None


class HologicSDKProvider:
    name = "hologic_sdk"
    experimental = False

    def __init__(self, sdk_path: Path, config: HologicSDKConfig) -> None:
        self._sdk_path = sdk_path
        self._config = config
        log.info("hologic_sdk_provider_init")

    def can_read(self, source: Path) -> bool:
        return False

    def load(self, source: Path) -> "QuantitativeMapBundle":
        raise NotImplementedError(
            "The Hologic SDK or its documentation has not been provided. "
            "To enable this provider: (1) obtain the Hologic SDK, "
            "(2) set HOLOGIC_SDK_PATH, (3) implement the binding in "
            "HologicSDKProvider.load() following the SDK API. "
            "See docs/dicom_hologic.md for what functions are needed."
        )

    def describe_capabilities(self) -> dict[str, object]:
        return {
            "status": "not_implemented",
            "blocked_by": "Hologic SDK not provided",
            "required_inputs": [
                "Hologic SDK library (DLL or shared object)",
                "Valid SDK license",
                "API documentation for reading P/R files and returning calibrated maps",
            ],
            "experimental": False,
        }
