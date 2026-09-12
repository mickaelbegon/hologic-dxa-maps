"""Structured logging configuration.

Uses structlog with JSON output in CI/production and human-readable output
in interactive sessions. PHI scrubbing is enforced at the processor level.
"""

from __future__ import annotations

import logging
import sys
from typing import Any

import structlog


_PHI_KEYS = frozenset({
    "PatientName", "PatientID", "PatientBirthDate", "PatientSex",
    "PatientAge", "PatientWeight", "PatientSize",
    "OtherPatientIDs", "OtherPatientNames",
    "InstitutionName", "InstitutionAddress",
    "ReferringPhysicianName", "PerformingPhysicianName",
    "RequestingPhysician", "OperatorsName",
    "AccessionNumber",
    # research pseudonym is OK; raw patient identifiers are not
})


def _scrub_phi(
    logger: logging.Logger,  # noqa: ARG001
    method: str,  # noqa: ARG001
    event_dict: dict[str, Any],
) -> dict[str, Any]:
    """Remove any PHI keys that accidentally entered a log event."""
    for key in _PHI_KEYS:
        if key in event_dict:
            event_dict[key] = "[REDACTED]"
    return event_dict


def configure_logging(level: str = "INFO", *, json: bool = False) -> None:
    """Call once at startup from the CLI entry point."""
    shared_processors: list[Any] = [
        structlog.contextvars.merge_contextvars,
        structlog.stdlib.add_log_level,
        structlog.stdlib.add_logger_name,
        structlog.processors.TimeStamper(fmt="iso"),
        _scrub_phi,
    ]

    if json or not sys.stderr.isatty():
        renderer: Any = structlog.processors.JSONRenderer()
    else:
        renderer = structlog.dev.ConsoleRenderer(colors=True)

    structlog.configure(
        processors=[
            *shared_processors,
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )

    formatter = structlog.stdlib.ProcessorFormatter(
        processor=renderer,
        foreign_pre_chain=shared_processors,
    )

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(getattr(logging, level.upper(), logging.INFO))


def get_logger(name: str) -> structlog.stdlib.BoundLogger:
    return structlog.get_logger(name)
