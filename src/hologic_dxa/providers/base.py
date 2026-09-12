"""Abstract provider protocol and module-level singleton registry."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import structlog

if TYPE_CHECKING:
    from hologic_dxa.maps.bundle import QuantitativeMapBundle

log = structlog.get_logger(__name__)


class QuantitativeDataUnavailableError(Exception):
    """Raised when required calibration, geometry, or data is absent.

    The message must explain exactly what is missing and what is needed to
    unblock the computation — a bare "data not found" is not acceptable.
    """


@runtime_checkable
class QuantitativeMapProvider(Protocol):
    """Interface all map providers must satisfy."""

    name: str
    experimental: bool

    def can_read(self, source: Path) -> bool:
        """Return True if this provider can handle the given source path."""
        ...

    def load(self, source: Path) -> "QuantitativeMapBundle":
        """Load and return a validated QuantitativeMapBundle.

        Raises
        ------
        QuantitativeDataUnavailableError
            When required calibration, geometry, or raw data is absent.
        """
        ...

    def describe_capabilities(self) -> dict[str, object]:
        """Return a dict describing what this provider can and cannot do."""
        ...


class ProviderRegistry:
    """Ordered registry of known providers; first match wins."""

    def __init__(self) -> None:
        self._providers: list[QuantitativeMapProvider] = []

    def register(self, provider: QuantitativeMapProvider) -> None:
        log.info(
            "provider_registered",
            name=provider.name,
            experimental=provider.experimental,
        )
        self._providers.append(provider)

    def find_provider(self, source: Path) -> QuantitativeMapProvider | None:
        for provider in self._providers:
            if provider.can_read(source):
                log.debug("provider_matched", name=provider.name)
                return provider
        log.debug("no_provider_matched", source_exists=source.exists())
        return None

    def list_capabilities(self) -> list[dict[str, object]]:
        return [
            {"name": p.name, "experimental": p.experimental, **p.describe_capabilities()}
            for p in self._providers
        ]


registry: ProviderRegistry = ProviderRegistry()
