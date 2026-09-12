"""Example: validate mass conservation of quantitative maps against SR results.

Requires:
- An HDF5 map bundle (from export-hdf5 or imported arrays)
- SR long-format CSV (from parse-sr + normalize-sr)

Usage:
    python examples/validate_maps.py output/maps/bundle.h5 output/sr/sr_long_format.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

from hologic_dxa.logging import configure_logging


def main() -> None:
    configure_logging("INFO")
    print(
        "This example requires a validated QuantitativeMapBundle (HDF5) and "
        "a normalized SR CSV. See docs/validation_protocol.md for the full workflow."
    )
    print("\nTo run validation:")
    print("  hologic-dxa validate-maps output/maps/bundle.h5 --sr output/sr/sr_long_format.csv")


if __name__ == "__main__":
    main()
