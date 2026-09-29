"""Whole-body attenuation images from the Hologic R file (tag 0023,1004).

EXPERIMENTAL: the layout and the per-plane pedestals were obtained by reverse
engineering (Horizon W, APEX 13.6; see docs/dicom_hologic.md section 8) and are
not validated for other devices or software versions.

Image convention (after parsing): axis 0 = lateral (row 0 = patient right),
axis 1 = head-to-foot (col 0 = top of head).
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np
from pydicom.dataset import Dataset

R_FILE_TAG = (0x0023, 0x1004)
IMAGE_RECORD_TYPES = (202, 224)
DESCRIPTOR_RECORD_TYPE = 221

# Empirical dark/pedestal level of each of the 6 energy planes.
PEDESTALS: tuple[float, ...] = (223.0, 372.0, 361.0, 622.0, 308.0, 474.0)
# Empirical native pixel pitch of the whole-body grid (lateral x head-to-foot).
DEFAULT_ROW_SIZE_M = 0.0057
DEFAULT_COL_SIZE_M = 0.0131


@dataclass(frozen=True)
class WholeBodyMu:
    """High- and low-energy attenuation maps, shape (lateral, head-to-foot)."""

    mu_h: np.ndarray
    mu_l: np.ndarray
    row_size_m: float = DEFAULT_ROW_SIZE_M
    col_size_m: float = DEFAULT_COL_SIZE_M


def parse_r_file_images(data: bytes) -> list[np.ndarray]:
    """Parse the R-file TLV stream into raw uint16 arrays of shape (planes, rows, cols)."""
    n, pos = len(data), 0
    records: list[tuple[int, bytes]] = []
    while pos + 6 <= n:
        rtype = struct.unpack_from("<H", data, pos)[0]
        length = struct.unpack_from("<I", data, pos + 2)[0]
        if length < 6 or pos + length > n:
            break
        records.append((rtype, data[pos + 6 : pos + length]))
        if rtype == 0:
            break
        pos += length

    images: list[np.ndarray] = []
    for i, (rtype, payload) in enumerate(records):
        if rtype != DESCRIPTOR_RECORD_TYPE or len(payload) < 26:
            continue
        n_rows = struct.unpack_from("<H", payload, 20)[0]
        n_planes = struct.unpack_from("<H", payload, 24)[0]
        if n_rows == 0 or n_planes == 0:
            continue
        for next_type, next_payload in records[i + 1 :]:
            if next_type not in IMAGE_RECORD_TYPES:
                break
            n_cols = len(next_payload) // (n_rows * n_planes * 2)
            if n_cols == 0:
                break
            raw = np.frombuffer(next_payload[: n_rows * n_planes * n_cols * 2], dtype="<u2")
            images.append(raw.reshape(n_cols, n_rows, n_planes).transpose(2, 1, 0))
    return images


def load_whole_body_mu(ds: Dataset) -> WholeBodyMu:
    """Return pedestal-corrected mu_H / mu_L from a Hologic whole-body dataset.

    Raises
    ------
    ValueError
        If tag (0023,1004) is absent or holds no parsable body image.
    """
    if R_FILE_TAG not in ds:
        raise ValueError(
            "Tag (0023,1004) (Hologic R file) is absent: cannot build whole-body "
            "attenuation maps. Is this a Hologic whole-body DXA DICOM?"
        )
    images = parse_r_file_images(bytes(ds[R_FILE_TAG].value))
    if not images:
        raise ValueError(
            "No body image found in R file (expected record 221 followed by 202/224)."
        )
    planes = np.stack([img.astype(np.float64) for img in images]).mean(axis=0)
    if planes.shape[0] < len(PEDESTALS):
        raise ValueError(f"Expected {len(PEDESTALS)} energy planes, found {planes.shape[0]}.")
    for p, pedestal in enumerate(PEDESTALS):
        planes[p] -= pedestal
    mu_h = np.maximum(planes[[0, 2, 4]].mean(axis=0), 0.0)
    mu_l = np.maximum(planes[[1, 3, 5]].mean(axis=0), 0.0)
    return WholeBodyMu(mu_h=mu_h, mu_l=mu_l)
