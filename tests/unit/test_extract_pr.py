"""Unit tests for hologic_dxa.dicom.extract_pr."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from hologic_dxa.dicom.extract_pr import (
    describe_pr,
    extract_pr_from_dataset,
    extract_pr_from_path,
)
from synthetic_data.generators import (
    make_hologic_archive,
    make_secondary_capture,
    save_to_tmpdir,
)

_RESEARCH_ID = "SUBJ001"
_FAKE_SHA = "a" * 64


# ---------------------------------------------------------------------------
# extract_pr_from_dataset
# ---------------------------------------------------------------------------

class TestExtractPrFromDataset:
    def test_happy_path_both_p_and_r(self, hologic_archive_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )

        assert manifest["source_dicom_sha256"] == _FAKE_SHA
        assert manifest["sop_instance_uid"]
        assert manifest["p_file"] is not None
        assert manifest["r_file"] is not None

        for key in ("p_file", "r_file"):
            entry = manifest[key]
            assert "path" in entry
            assert "sha256" in entry
            assert "extracted_length" in entry
            assert Path(entry["path"]).exists()

        p_bytes = Path(manifest["p_file"]["path"]).read_bytes()
        assert hashlib.sha256(p_bytes).hexdigest() == manifest["p_file"]["sha256"]

    def test_p_only(self, hologic_archive_p_only_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_p_only_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        assert manifest["p_file"] is not None
        assert manifest["r_file"] is None

    def test_empty_blob_writes_file_and_warns(self, hologic_archive_empty_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_empty_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        p = manifest["p_file"]
        assert p is not None
        assert p["extracted_length"] == 0
        assert Path(p["path"]).exists()
        assert any("empty" in w.lower() for w in manifest["warnings"])

    def test_padded_blob_warns_no_exception(self, hologic_archive_padded_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_padded_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        assert manifest["p_file"] is not None
        warns = manifest["warnings"]
        assert any("1 byte" in w or "padding" in w.lower() for w in warns)

    def test_truncated_raises_value_error(self, hologic_archive_truncated_ds, tmp_path):
        with pytest.raises(ValueError, match="mismatch"):
            extract_pr_from_dataset(
                hologic_archive_truncated_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
            )

    def test_no_pr_tags_raises_runtime_error(self, tmp_path):
        ds = make_secondary_capture()
        with pytest.raises(RuntimeError, match="neither P nor R"):
            extract_pr_from_dataset(ds, _FAKE_SHA, tmp_path, _RESEARCH_ID)

    def test_output_dir_created_if_absent(self, hologic_archive_ds, tmp_path):
        nested = tmp_path / "deep" / "nested"
        assert not nested.exists()
        extract_pr_from_dataset(hologic_archive_ds, _FAKE_SHA, nested, _RESEARCH_ID)
        assert nested.is_dir()

    def test_filename_starts_with_research_id(self, hologic_archive_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        p_name = Path(manifest["p_file"]["path"]).name
        assert p_name.startswith(_RESEARCH_ID)

    def test_declared_length_recorded(self, hologic_archive_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        assert manifest["p_file"]["declared_length"] is not None
        assert manifest["p_file"]["declared_length"] == len(b"FAKE_P_FILE_CONTENT")

    def test_warnings_list_present(self, hologic_archive_ds, tmp_path):
        manifest = extract_pr_from_dataset(
            hologic_archive_ds, _FAKE_SHA, tmp_path, _RESEARCH_ID
        )
        assert isinstance(manifest["warnings"], list)

    def test_custom_p_content_written_verbatim(self, tmp_path):
        p_bytes = b"\x00\x01\x02\x03" * 16
        ds = make_hologic_archive(p_data=p_bytes)
        manifest = extract_pr_from_dataset(ds, _FAKE_SHA, tmp_path, _RESEARCH_ID)
        written = Path(manifest["p_file"]["path"]).read_bytes()
        assert written == p_bytes


# ---------------------------------------------------------------------------
# extract_pr_from_path
# ---------------------------------------------------------------------------

class TestExtractPrFromPath:
    def test_writes_manifest_json(self, tmp_path):
        ds = make_hologic_archive()
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        dicom_path = save_to_tmpdir(ds, input_dir, "archive.dcm")
        output_dir = tmp_path / "output"

        manifest = extract_pr_from_path(dicom_path, output_dir, _RESEARCH_ID)

        assert manifest is not None
        manifests = list(output_dir.glob("*_manifest.json"))
        assert len(manifests) == 1

        written = json.loads(manifests[0].read_text("utf-8"))
        assert written["sop_instance_uid"] == manifest["sop_instance_uid"]
        assert written["p_file"] is not None

    def test_returns_none_for_non_archive(self, tmp_path):
        ds = make_secondary_capture()
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        dicom_path = save_to_tmpdir(ds, input_dir, "sc.dcm")
        result = extract_pr_from_path(dicom_path, tmp_path / "output", _RESEARCH_ID)
        assert result is None

    def test_source_sha256_in_manifest(self, tmp_path):
        ds = make_hologic_archive()
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        dicom_path = save_to_tmpdir(ds, input_dir, "archive.dcm")
        output_dir = tmp_path / "output"

        manifest = extract_pr_from_path(dicom_path, output_dir, _RESEARCH_ID)
        assert manifest is not None

        expected_sha = hashlib.sha256(dicom_path.read_bytes()).hexdigest()
        assert manifest["source_dicom_sha256"] == expected_sha


# ---------------------------------------------------------------------------
# describe_pr
# ---------------------------------------------------------------------------

class TestDescribePr:
    def _write_manifest(
        self,
        tmp_path: Path,
        p_bytes: bytes,
        r_bytes: bytes | None = None,
        corrupt_p_sha: bool = False,
    ) -> Path:
        p_path = tmp_path / "p.bin"
        p_path.write_bytes(p_bytes)
        p_sha = hashlib.sha256(p_bytes).hexdigest()
        if corrupt_p_sha:
            p_sha = "badhash"

        r_entry = None
        if r_bytes is not None:
            r_path = tmp_path / "r.bin"
            r_path.write_bytes(r_bytes)
            r_entry = {
                "path": str(r_path),
                "declared_length": len(r_bytes),
                "extracted_length": len(r_bytes),
                "sha256": hashlib.sha256(r_bytes).hexdigest(),
            }

        manifest = {
            "source_dicom_sha256": "deadbeef",
            "sop_instance_uid": "1.2.3.4.5",
            "encoding_scheme": "HOLOGIC_ENCODING_V1",
            "p_file": {
                "path": str(p_path),
                "declared_length": len(p_bytes),
                "extracted_length": len(p_bytes),
                "sha256": p_sha,
            },
            "r_file": r_entry,
            "warnings": [],
        }
        mp = tmp_path / "manifest.json"
        mp.write_text(json.dumps(manifest), encoding="utf-8")
        return mp

    def test_returns_forensic_fields(self, tmp_path):
        mp = self._write_manifest(tmp_path, b"FAKE_P" * 100, b"FAKE_R" * 100)
        desc = describe_pr(mp)

        assert desc["sop_instance_uid"] == "1.2.3.4.5"
        assert desc["encoding_scheme"] == "HOLOGIC_ENCODING_V1"

        for key in ("p_file", "r_file"):
            entry = desc[key]
            assert entry is not None
            assert isinstance(entry["size_bytes"], int)
            assert isinstance(entry["byte_entropy_bits"], float)
            assert isinstance(entry["first_16_bytes_hex"], str)
            assert isinstance(entry["printable_ascii_ratio"], float)
            assert entry["sha256_verified"] is True

    def test_sha256_mismatch_flagged(self, tmp_path):
        mp = self._write_manifest(tmp_path, b"CONTENT", b"R", corrupt_p_sha=True)
        desc = describe_pr(mp)
        assert desc["p_file"]["sha256_verified"] is False
        assert desc["r_file"]["sha256_verified"] is True

    def test_r_file_none_in_manifest(self, tmp_path):
        mp = self._write_manifest(tmp_path, b"P_ONLY", r_bytes=None)
        desc = describe_pr(mp)
        assert desc["p_file"] is not None
        assert desc["r_file"] is None

    def test_missing_manifest_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            describe_pr(tmp_path / "nonexistent_manifest.json")

    def test_missing_binary_raises_file_not_found(self, tmp_path):
        manifest = {
            "source_dicom_sha256": "x",
            "sop_instance_uid": "1.2.3",
            "encoding_scheme": "",
            "p_file": {
                "path": str(tmp_path / "missing.bin"),
                "declared_length": 10,
                "extracted_length": 10,
                "sha256": "abc",
            },
            "r_file": None,
            "warnings": [],
        }
        mp = tmp_path / "manifest.json"
        mp.write_text(json.dumps(manifest), encoding="utf-8")
        with pytest.raises(FileNotFoundError):
            describe_pr(mp)

    def test_first_16_bytes_hex_length(self, tmp_path):
        mp = self._write_manifest(tmp_path, b"\xde\xad\xbe\xef" * 4)
        desc = describe_pr(mp)
        assert len(desc["p_file"]["first_16_bytes_hex"]) == 32  # 16 bytes = 32 hex chars

    def test_note_field_present(self, tmp_path):
        mp = self._write_manifest(tmp_path, b"ANY_CONTENT")
        desc = describe_pr(mp)
        assert "note" in desc["p_file"]
