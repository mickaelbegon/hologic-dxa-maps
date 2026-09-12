"""Extraction of Hologic P/R binary files from DICOM archives.

Security constraints (non-negotiable)
--------------------------------------
* The filename stored in tag (0023,1001) is NEVER used as a filesystem path.
  It is treated as an opaque string for informational logging only.
* Extracted bytes are written verbatim; no transformation is applied.
* Output filenames are derived exclusively from the research ID and the SOP
  Instance UID, not from any tag value that arrives in the DICOM file.

Format agnosticism
------------------
``describe_pr`` performs forensic inspection (entropy, magic bytes, ASCII
ratio) but deliberately avoids claiming to know the binary format unless
concrete evidence is present. Hologic P/R format documentation is not
publicly available; any format claim must be backed by external references.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pydicom
import structlog

from hologic_dxa.dicom.private_tags import HologicTag, inspect_private_tags
from hologic_dxa.provenance import sha256_bytes

logger = structlog.get_logger(__name__)

# Known magic-byte signatures.  Extend as evidence accumulates.
_MAGIC_SIGNATURES: list[tuple[bytes, str]] = [
    (b"PK\x03\x04", "ZIP archive"),
    (b"\x1f\x8b", "gzip"),
    (b"BZh", "bzip2"),
    (b"\xfd7zXZ\x00", "XZ/LZMA"),
    (b"RIFF", "RIFF (WAV/AVI)"),
    (b"\x89PNG\r\n\x1a\n", "PNG image"),
    (b"\xff\xd8\xff", "JPEG image"),
    (b"%PDF", "PDF document"),
    (b"MZ", "Windows PE executable"),
    (b"\x7fELF", "ELF executable"),
    (b"OggS", "Ogg container"),
    (b"fLaC", "FLAC audio"),
    (b"\x00\x00\x00\x0cftyp", "MPEG-4 container"),
]

_UID_SAFE_RE = re.compile(r"[^A-Za-z0-9]")


def _sanitize_uid(uid: str, max_len: int = 64) -> str:
    """Replace non-alphanumeric characters with underscores, truncate to *max_len*."""
    return _UID_SAFE_RE.sub("_", uid)[:max_len]


def _sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _byte_entropy(data: bytes) -> float:
    """Shannon entropy of the byte-value distribution, in bits per byte.

    Uses ``scipy.stats.entropy`` on the 256-bin byte histogram.
    Returns 0.0 for empty data.
    """
    if not data:
        return 0.0
    from scipy.stats import entropy as scipy_entropy  # deferred to avoid hard dep at import

    counts = [0] * 256
    for b in data:
        counts[b] += 1
    total = len(data)
    probs = [c / total for c in counts if c > 0]
    # scipy_entropy uses natural log by default; request log base 2
    return float(scipy_entropy(probs, base=2))


def _detect_magic(data: bytes) -> str | None:
    """Return the first matching magic-byte description, or None."""
    for magic, description in _MAGIC_SIGNATURES:
        if data[: len(magic)] == magic:
            return description
    return None


def _printable_ascii_ratio(data: bytes) -> float:
    """Fraction of bytes that are printable ASCII (0x20–0x7E)."""
    if not data:
        return 0.0
    count = sum(0x20 <= b <= 0x7E for b in data)
    return count / len(data)


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def extract_pr_from_dataset(
    ds: pydicom.Dataset,
    source_sha256: str,
    output_dir: Path,
    research_id: str,
) -> dict[str, Any]:
    """Extract P and R binary files from a pydicom Dataset.

    Parameters
    ----------
    ds:
        Dataset already loaded by pydicom (pixels must be present).
    source_sha256:
        Pre-computed SHA-256 of the source DICOM file (caller's responsibility).
    output_dir:
        Directory where the extracted files will be written.  Created if absent.
    research_id:
        Pseudonymised subject identifier (see :class:`PseudonymRegistry`).

    Returns
    -------
    dict
        Manifest as specified by the module-level schema.

    Raises
    ------
    RuntimeError
        If neither P nor R data is present in the dataset.
    ValueError
        If extracted length deviates from the declared length by more than 1 byte.
    """
    inspection = inspect_private_tags(ds)
    sop_uid = inspection.sop_instance_uid
    uid_safe = _sanitize_uid(sop_uid)

    private_creator = inspection.private_creator or ""
    encoding_scheme = inspection.encoding_scheme or ""

    output_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = list(inspection.warnings)
    manifest: dict[str, Any] = {
        "source_dicom_sha256": source_sha256,
        "sop_instance_uid": sop_uid,
        "encoding_scheme": encoding_scheme,
        "private_creator": private_creator,
        "p_file": None,
        "r_file": None,
        "warnings": warnings,
    }

    p_tag = HologicTag.P_FILE_DATA.value
    r_tag = HologicTag.R_FILE_DATA.value
    p_len_tag = HologicTag.P_FILE_LENGTH.value
    r_len_tag = HologicTag.R_FILE_LENGTH.value

    has_p = p_tag in ds
    has_r = r_tag in ds

    if not has_p and not has_r:
        raise RuntimeError(
            f"Dataset {sop_uid} contains neither P nor R file data tags. "
            "Call inspect_private_tags first to confirm PR_PRESENT before extracting."
        )

    def _extract_one(
        data_tag: object,
        len_tag: object,
        suffix: str,
    ) -> dict[str, Any] | None:
        if data_tag not in ds:
            return None

        raw = ds[data_tag].value  # type: ignore[index]
        if not isinstance(raw, (bytes, bytearray)):
            warnings.append(
                f"{suffix} file data tag value is not bytes (got {type(raw).__name__}). "
                "Skipping."
            )
            return None

        blob: bytes = bytes(raw)
        extracted_len = len(blob)

        if extracted_len == 0:
            warnings.append(f"{suffix} file blob is empty (0 bytes).")

        declared_len: int | None = None
        if len_tag in ds:
            try:
                declared_len = int(ds[len_tag].value)  # type: ignore[index]
            except (ValueError, TypeError) as exc:
                warnings.append(f"Could not parse declared {suffix} length: {exc}")

        if declared_len is not None:
            diff = abs(extracted_len - declared_len)
            if diff > 1:
                raise ValueError(
                    f"{suffix} file length mismatch: declared={declared_len}, "
                    f"extracted={extracted_len} (difference={diff}). "
                    "Allowing ±1 byte for DICOM odd-length padding only."
                )
            if diff == 1:
                warnings.append(
                    f"{suffix} file length differs by 1 byte from declared value "
                    f"(declared={declared_len}, extracted={extracted_len}). "
                    "Consistent with DICOM odd-length padding."
                )

        filename = f"{research_id}_{uid_safe}_{suffix}.bin"
        out_path = output_dir / filename
        out_path.write_bytes(blob)

        digest = sha256_bytes(blob)
        logger.info(
            "pr_extraction.wrote",
            suffix=suffix,
            filename=filename,
            extracted_length=extracted_len,
            sha256=digest,
        )

        return {
            "path": str(out_path),
            "declared_length": declared_len,
            "extracted_length": extracted_len,
            "sha256": digest,
        }

    manifest["p_file"] = _extract_one(p_tag, p_len_tag, "P")
    manifest["r_file"] = _extract_one(r_tag, r_len_tag, "R")

    return manifest


def extract_pr_from_path(
    dicom_path: Path,
    output_dir: Path,
    research_id: str,
) -> dict[str, Any] | None:
    """Load a DICOM file and extract its P/R data.

    Parameters
    ----------
    dicom_path:
        Path to the DICOM file.  The path is only used for I/O; the filename
        is NEVER used to construct output paths or infer format.
    output_dir:
        Directory for extracted binaries and the manifest JSON.
    research_id:
        Pseudonymised subject identifier.

    Returns
    -------
    dict | None
        Manifest dict, or ``None`` if the file contains no P/R tags.
    """
    logger.info("pr_extraction.loading", path=str(dicom_path))

    source_sha256 = _sha256_of_file(dicom_path)
    ds = pydicom.dcmread(str(dicom_path))

    inspection = inspect_private_tags(ds)
    from hologic_dxa.dicom.private_tags import PRDiagnosis

    if inspection.diagnosis not in (PRDiagnosis.PR_PRESENT, PRDiagnosis.PR_PARTIAL):
        logger.debug(
            "pr_extraction.no_pr_tags",
            diagnosis=inspection.diagnosis.value,
        )
        return None

    manifest = extract_pr_from_dataset(ds, source_sha256, output_dir, research_id)

    # Write the manifest alongside the extracted files.
    sop_uid_safe = _sanitize_uid(inspection.sop_instance_uid)
    manifest_path = output_dir / f"{research_id}_{sop_uid_safe}_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, default=str),
        encoding="utf-8",
    )
    logger.info("pr_extraction.manifest_written", path=str(manifest_path))

    return manifest


# ---------------------------------------------------------------------------
# Forensic description
# ---------------------------------------------------------------------------

def describe_pr(manifest_path: Path) -> dict[str, Any]:
    """Return a forensic description of the P/R files named in a manifest.

    Measurements include byte count, Shannon entropy, first-16-bytes hex dump,
    detected magic bytes, and printable-ASCII ratio. No interpretation of the
    binary format is made beyond what the magic-byte table can prove.

    Parameters
    ----------
    manifest_path:
        Path to a JSON manifest produced by :func:`extract_pr_from_path`.

    Returns
    -------
    dict
        Forensic summary keyed by ``"p_file"`` and ``"r_file"``.

    Raises
    ------
    FileNotFoundError
        If the manifest or one of the referenced binary files does not exist.
    ValueError
        If the manifest is malformed.
    """
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))

    description: dict[str, Any] = {
        "manifest_path": str(manifest_path),
        "sop_instance_uid": raw.get("sop_instance_uid", ""),
        "encoding_scheme": raw.get("encoding_scheme", ""),
        "p_file": None,
        "r_file": None,
    }

    for key in ("p_file", "r_file"):
        entry = raw.get(key)
        if entry is None:
            continue

        bin_path = Path(entry["path"])
        if not bin_path.exists():
            raise FileNotFoundError(
                f"Binary file referenced by manifest not found: {bin_path}"
            )

        data = bin_path.read_bytes()
        size = len(data)
        first16 = data[:16].hex()
        entropy = _byte_entropy(data)
        magic = _detect_magic(data)
        ascii_ratio = _printable_ascii_ratio(data)

        # Re-verify SHA-256 to catch corruption after extraction.
        actual_sha256 = sha256_bytes(data)
        integrity_ok = actual_sha256 == entry.get("sha256", "")

        suffix_desc: dict[str, Any] = {
            "path": str(bin_path),
            "size_bytes": size,
            "sha256_verified": integrity_ok,
            "byte_entropy_bits": round(entropy, 4),
            "first_16_bytes_hex": first16,
            "detected_magic": magic,
            "printable_ascii_ratio": round(ascii_ratio, 4),
            "note": (
                "Format is UNKNOWN unless detected_magic provides proof. "
                "Hologic P/R format documentation is not publicly available."
            ),
        }
        description[key] = suffix_desc

    return description
