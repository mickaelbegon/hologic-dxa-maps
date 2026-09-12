"""Example: inspect a directory of DICOM files and print the inventory summary.

Usage:
    python examples/inspect_dicom.py /path/to/dicoms

Or use the CLI:
    hologic-dxa inspect /path/to/dicoms --output output/inventory
"""

from __future__ import annotations

import sys
from pathlib import Path

from hologic_dxa.dicom.inventory import build_inventory, write_inventory
from hologic_dxa.logging import configure_logging


def main() -> None:
    configure_logging("INFO")

    if len(sys.argv) < 2:
        print("Usage: python inspect_dicom.py <dicom_directory>")
        sys.exit(1)

    root = Path(sys.argv[1])
    if not root.is_dir():
        print(f"Error: '{root}' is not a directory.")
        sys.exit(1)

    records = build_inventory(root)
    output_dir = Path("output/inventory")
    write_inventory(records, output_dir)

    dicom_count = sum(r.is_dicom for r in records)
    print(f"\nProcessed {len(records)} files, {dicom_count} valid DICOM.")
    print(f"Results written to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
