"""Provenance tracking and pseudonymisation registry for the DXA pipeline.

Provenance records capture everything needed to reproduce or audit a pipeline
run without storing any patient-identifying information.
"""

from __future__ import annotations

import hashlib
import json
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import structlog

from hologic_dxa import __version__

logger = structlog.get_logger(__name__)


@dataclass
class Provenance:
    """Immutable record of a single pipeline operation.

    All fields are set at construction time; timestamps are stored as UTC.
    """

    pipeline_version: str
    timestamp_utc: datetime
    input_hashes: dict[str, str]
    operations: list[str]
    source_sop_instance_uid: str
    study_instance_uid: str
    series_instance_uid: str
    device_manufacturer: str
    device_model: str
    apex_version: str | None
    transfer_syntax_uid: str

    @classmethod
    def create(
        cls,
        *,
        input_hashes: dict[str, str],
        source_sop_instance_uid: str,
        study_instance_uid: str,
        series_instance_uid: str,
        device_manufacturer: str,
        device_model: str,
        apex_version: str | None,
        transfer_syntax_uid: str,
        operations: list[str] | None = None,
    ) -> "Provenance":
        """Construct a Provenance stamped with the current UTC time."""
        return cls(
            pipeline_version=__version__,
            timestamp_utc=datetime.now(tz=timezone.utc),
            input_hashes=dict(input_hashes),
            operations=list(operations or []),
            source_sop_instance_uid=source_sop_instance_uid,
            study_instance_uid=study_instance_uid,
            series_instance_uid=series_instance_uid,
            device_manufacturer=device_manufacturer,
            device_model=device_model,
            apex_version=apex_version,
            transfer_syntax_uid=transfer_syntax_uid,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable representation of this record."""
        return {
            "pipeline_version": self.pipeline_version,
            "timestamp_utc": self.timestamp_utc.isoformat(),
            "input_hashes": dict(self.input_hashes),
            "operations": list(self.operations),
            "source_sop_instance_uid": self.source_sop_instance_uid,
            "study_instance_uid": self.study_instance_uid,
            "series_instance_uid": self.series_instance_uid,
            "device_manufacturer": self.device_manufacturer,
            "device_model": self.device_model,
            "apex_version": self.apex_version,
            "transfer_syntax_uid": self.transfer_syntax_uid,
        }

    def record_operation(self, description: str) -> None:
        """Append an operation description to the ordered log."""
        self.operations.append(description)


def sha256_file(path: Path) -> str:
    """Return the hex SHA-256 digest of *path* without loading it all into RAM."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    """Return the hex SHA-256 digest of *data*."""
    return hashlib.sha256(data).hexdigest()


class PseudonymRegistry:
    """Bijective map between patient IDs and opaque research IDs.

    The registry is backed by a JSON file so research IDs persist across runs.
    Patient IDs are NEVER logged; only research IDs appear in log events.

    Research IDs have the form ``SUB_{8-char hex prefix}``, e.g. ``SUB_A3F1C2B0``.
    """

    _ID_PREFIX = "SUB_"

    def __init__(self, registry_path: Path) -> None:
        self._path = registry_path
        # _forward: patient_id → research_id  (never logged)
        self._forward: dict[str, str] = {}
        # _reverse: research_id → True  (existence check only)
        self._reverse: set[str] = set()
        if registry_path.exists():
            self.load()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_or_create(self, patient_id: str) -> str:
        """Return an existing research ID or allocate a new one.

        The *patient_id* argument is treated as a secret and is never
        forwarded to any logger or external system.
        """
        if patient_id in self._forward:
            return self._forward[patient_id]

        research_id = self._generate_unique_id()
        self._forward[patient_id] = research_id
        self._reverse.add(research_id)
        self.save()
        logger.info("pseudonym.allocated", research_id=research_id)
        return research_id

    def load(self) -> None:
        """Reload the registry from disk, replacing in-memory state."""
        raw: dict[str, str] = json.loads(self._path.read_text(encoding="utf-8"))

        # Validate bijectivity before replacing state
        seen_research: set[str] = set()
        for research_id in raw.values():
            if research_id in seen_research:
                raise ValueError(
                    f"Registry at '{self._path}' contains a duplicate research_id "
                    f"'{research_id}'. The file may have been manually edited or merged "
                    "incorrectly. Fix the collision before proceeding."
                )
            seen_research.add(research_id)

        self._forward = dict(raw)
        self._reverse = seen_research
        logger.debug("pseudonym_registry.loaded", entry_count=len(self._forward))

    def save(self) -> None:
        """Persist the current registry to disk atomically."""
        self._path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self._path.with_suffix(".tmp")
        tmp.write_text(
            json.dumps(self._forward, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        tmp.replace(self._path)
        logger.debug("pseudonym_registry.saved", entry_count=len(self._forward))

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _generate_unique_id(self) -> str:
        for _ in range(1_000):
            candidate = self._ID_PREFIX + uuid.uuid4().hex[:8].upper()
            if candidate not in self._reverse:
                return candidate
        raise RuntimeError(
            "Failed to generate a unique research ID after 1000 attempts. "
            "This should never happen with a UUID4 source."
        )
