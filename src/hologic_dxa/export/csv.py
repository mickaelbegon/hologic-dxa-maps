"""CSV export for mass properties and validation results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def export_mass_properties_csv(
    mass_properties: dict[str, dict],
    output_path: Path,
) -> None:
    """Write mass properties (mass, CoM, inertia per region) to CSV.

    Parameters
    ----------
    mass_properties:
        Output of maps.integration.compute_mass_properties().
    output_path:
        Destination CSV file.
    """
    rows = []
    for region, props in mass_properties.items():
        row = {"region": region}
        row.update(props)
        rows.append(row)

    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def export_sr_long_format(
    measurements: list[dict],
    output_path: Path,
) -> None:
    """Write SR measurements to long-format CSV."""
    df = pd.DataFrame(measurements)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
