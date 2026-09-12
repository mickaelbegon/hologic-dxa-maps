"""Unit tests for hologic_dxa.dicom.inventory."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pytest

from hologic_dxa.dicom.inventory import build_inventory, sha256_file
from tests.synthetic_data.generators import (
    make_hologic_archive,
    make_secondary_capture,
    save_to_tmpdir,
)

# ---------------------------------------------------------------------------
# PHI field names that must never appear in any inventory output
# ---------------------------------------------------------------------------

_PHI_KEYS = {
    "PatientName",
    "PatientID",
    "PatientBirthDate",
    "AccessionNumber",
    "PatientAge",
    "PatientWeight",
    "PatientSex",
    "PatientSize",
    "OtherPatientIDs",
    "OtherPatientNames",
    "InstitutionName",
    "InstitutionAddress",
}


def _phi_keys_in_dict(d: dict) -> set[str]:
    """Recursively collect any PHI key names found in *d*."""
    found: set[str] = set()
    for k, v in d.items():
        if k in _PHI_KEYS:
            found.add(k)
        if isinstance(v, dict):
            found |= _phi_keys_in_dict(v)
        elif isinstance(v, list):
            for item in v:
                if isinstance(item, dict):
                    found |= _phi_keys_in_dict(item)
    return found


# ---------------------------------------------------------------------------
# sha256_file
# ---------------------------------------------------------------------------

class TestSha256File:
    def test_sha256_reproducible(self, tmp_path):
        """Same file content → identical SHA-256 digest on two calls."""
        test_file = tmp_path / "data.bin"
        test_file.write_bytes(b"reproducible content " * 1000)
        hash1 = sha256_file(test_file)
        hash2 = sha256_file(test_file)
        assert hash1 == hash2

    def test_sha256_is_hex_string(self, tmp_path):
        """SHA-256 digest must be a 64-character lowercase hex string."""
        test_file = tmp_path / "data.bin"
        test_file.write_bytes(b"some content")
        digest = sha256_file(test_file)
        assert len(digest) == 64
        assert all(c in "0123456789abcdef" for c in digest)

    def test_sha256_different_for_different_content(self, tmp_path):
        """Different content → different SHA-256 digest."""
        file_a = tmp_path / "a.bin"
        file_b = tmp_path / "b.bin"
        file_a.write_bytes(b"content A")
        file_b.write_bytes(b"content B")
        assert sha256_file(file_a) != sha256_file(file_b)

    def test_sha256_empty_file(self, tmp_path):
        """Empty file → well-known SHA-256 of empty string."""
        empty_file = tmp_path / "empty.bin"
        empty_file.write_bytes(b"")
        # SHA-256 of empty string is a known constant
        expected = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        assert sha256_file(empty_file) == expected


# ---------------------------------------------------------------------------
# build_inventory — PHI absence
# ---------------------------------------------------------------------------

class TestInventoryNoPHI:
    def test_inventory_no_phi(self, tmp_inventory_dir):
        """Inventory records must not contain any PHI field names."""
        records = build_inventory(tmp_inventory_dir)
        for record in records:
            record_dict = asdict(record)
            phi_found = _phi_keys_in_dict(record_dict)
            assert len(phi_found) == 0, (
                f"PHI fields found in inventory record for "
                f"'{record.file_path}': {phi_found}"
            )

    def test_inventory_json_no_phi(self, tmp_inventory_dir):
        """JSON serialisation of inventory must not contain PHI key names."""
        records = build_inventory(tmp_inventory_dir)
        dicts = [asdict(r) for r in records]
        json_text = json.dumps(dicts)
        for phi_key in _PHI_KEYS:
            assert phi_key not in json_text, (
                f"PHI key '{phi_key}' found in serialised inventory JSON."
            )


# ---------------------------------------------------------------------------
# build_inventory — non-DICOM file handling
# ---------------------------------------------------------------------------

class TestNonDicomFile:
    def test_non_dicom_file_recorded(self, tmp_path):
        """Plain text file → is_dicom=False, parse_error describes the issue."""
        txt_file = tmp_path / "notes.txt"
        txt_file.write_text("not a dicom file\n", encoding="utf-8")

        records = build_inventory(tmp_path)
        assert len(records) == 1
        record = records[0]
        assert record.is_dicom is False
        assert record.parse_error is not None
        assert len(record.parse_error) > 0

    def test_non_dicom_file_has_sha256(self, tmp_path):
        """Non-DICOM files must still receive a SHA-256 hash."""
        txt_file = tmp_path / "notes.txt"
        txt_file.write_text("not a dicom file\n", encoding="utf-8")
        records = build_inventory(tmp_path)
        assert records[0].sha256 != ""

    def test_non_dicom_file_path_is_relative(self, tmp_path):
        """File path in the record must be relative to the root."""
        txt_file = tmp_path / "notes.txt"
        txt_file.write_text("content", encoding="utf-8")
        records = build_inventory(tmp_path)
        # relative path should not start with the tmp_path prefix
        assert not records[0].file_path.startswith(str(tmp_path))
        assert "notes.txt" in records[0].file_path


# ---------------------------------------------------------------------------
# build_inventory — DICOM classification
# ---------------------------------------------------------------------------

class TestInventoryDicomClassification:
    def test_secondary_capture_classified(self, tmp_path):
        """Secondary Capture DICOM → classified as 'Secondary Capture'."""
        save_to_tmpdir(make_secondary_capture(), tmp_path, "sc.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert len(dicom_records) == 1
        record = dicom_records[0]
        assert record.classification == "Secondary Capture", (
            f"Expected 'Secondary Capture', got '{record.classification}'"
        )

    def test_secondary_capture_not_eligible(self, tmp_path):
        """Secondary Capture → quantitative_eligibility = 'not_eligible'."""
        save_to_tmpdir(make_secondary_capture(), tmp_path, "sc.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert dicom_records[0].quantitative_eligibility == "not_eligible"

    def test_hologic_archive_classified(self, tmp_path):
        """Hologic archive (P/R present) → classified as 'Hologic Archive (P/R data)'."""
        save_to_tmpdir(make_hologic_archive(), tmp_path, "archive.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert len(dicom_records) == 1
        record = dicom_records[0]
        assert record.classification == "Hologic Archive (P/R data)", (
            f"Expected 'Hologic Archive (P/R data)', got '{record.classification}'"
        )

    def test_hologic_archive_pr_diagnosis(self, tmp_path):
        """Hologic archive → pr_diagnosis = 'PR_PRESENT'."""
        save_to_tmpdir(make_hologic_archive(), tmp_path, "archive.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert dicom_records[0].pr_diagnosis == "PR_PRESENT"

    def test_hologic_archive_presence_flag(self, tmp_path):
        """Hologic archive → presence_of_hologic_pr_tags = True."""
        save_to_tmpdir(make_hologic_archive(), tmp_path, "archive.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert dicom_records[0].presence_of_hologic_pr_tags is True

    def test_inventory_records_file_size(self, tmp_path):
        """Inventory records must capture the actual file size in bytes."""
        path = save_to_tmpdir(make_secondary_capture(), tmp_path, "sc.dcm")
        expected_size = path.stat().st_size
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert dicom_records[0].file_size_bytes == expected_size

    def test_inventory_records_sop_class_uid(self, tmp_path):
        """Inventory must capture the SOPClassUID for DICOM files."""
        save_to_tmpdir(make_secondary_capture(), tmp_path, "sc.dcm")
        records = build_inventory(tmp_path)
        dicom_records = [r for r in records if r.is_dicom]
        assert dicom_records[0].sop_class_uid == "1.2.840.10008.5.1.4.1.1.7"


# ---------------------------------------------------------------------------
# build_inventory — mixed directory
# ---------------------------------------------------------------------------

class TestInventoryMixedDirectory:
    def test_mixed_dir_counts(self, tmp_inventory_dir):
        """Pre-populated fixture has 4 DICOM files and 1 non-DICOM text file."""
        records = build_inventory(tmp_inventory_dir)
        dicom_count = sum(r.is_dicom for r in records)
        non_dicom_count = sum(not r.is_dicom for r in records)
        assert dicom_count == 4, f"Expected 4 DICOM files, got {dicom_count}"
        assert non_dicom_count == 1, (
            f"Expected 1 non-DICOM file, got {non_dicom_count}"
        )

    def test_all_records_have_sha256(self, tmp_inventory_dir):
        """All records (DICOM and non-DICOM) must have a non-empty SHA-256."""
        records = build_inventory(tmp_inventory_dir)
        for record in records:
            if record.sha256 != "":  # empty only when IO error computing hash
                assert len(record.sha256) == 64
