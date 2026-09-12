"""Example: extract Hologic P and R files from archive DICOM objects.

Usage:
    python examples/extract_pr.py /path/to/dicoms

Or use the CLI:
    hologic-dxa extract-pr /path/to/dicoms --output output/extracted
"""

from __future__ import annotations

import sys
from pathlib import Path

from hologic_dxa.logging import configure_logging


def main() -> None:
    configure_logging("INFO")
    print(
        "P/R extraction extracts raw binary blobs from Hologic archive DICOMs.\n"
        "The extracted bytes are NOT decoded — decoding requires Hologic SDK or docs.\n"
        "See docs/dicom_hologic.md and docs/limitations.md for details.\n"
    )
    print("CLI command:")
    print("  hologic-dxa extract-pr /path/to/dicoms --output output/extracted")
    print("\nThis will:")
    print("  1. Detect DICOMs with Hologic private tags (group 0023)")
    print("  2. Extract P and R file bytes verbatim")
    print("  3. Write {research_id}_{uid}_P.bin and {research_id}_{uid}_R.bin")
    print("  4. Compute SHA-256 checksums")
    print("  5. Write a JSON manifest per object")


if __name__ == "__main__":
    main()
