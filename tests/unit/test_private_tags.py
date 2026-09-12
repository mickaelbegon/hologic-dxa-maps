"""Unit tests for hologic_dxa.dicom.private_tags.inspect_private_tags."""

from __future__ import annotations

import pytest

from hologic_dxa.dicom.private_tags import PRDiagnosis, inspect_private_tags
from synthetic_data.generators import (
    make_empty_pr_archive,
    make_hologic_archive,
    make_secondary_capture,
    make_truncated_pr_archive,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PHI_STRINGS = {
    "PatientName",
    "PatientID",
    "PatientBirthDate",
    "AccessionNumber",
    "PatientAge",
    "PatientWeight",
    "InstitutionName",
}


def _result_contains_phi(result: object) -> bool:
    """Return True if any PHI-like key appears anywhere in the result repr."""
    text = repr(result)
    return any(phi in text for phi in _PHI_STRINGS)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSecondaryCapture:
    def test_secondary_capture_returns_no_pr(self, secondary_capture_ds):
        """Secondary Capture: diagnosis is SECONDARY_CAPTURE_ONLY or ARCHIVE_DATA_ABSENT."""
        result = inspect_private_tags(secondary_capture_ds)
        assert result.diagnosis in (
            PRDiagnosis.SECONDARY_CAPTURE_ONLY,
            PRDiagnosis.ARCHIVE_DATA_ABSENT,
        ), f"Unexpected diagnosis: {result.diagnosis}"

    def test_secondary_capture_has_no_p_data(self, secondary_capture_ds):
        """Secondary Capture must not report any P or R file data."""
        result = inspect_private_tags(secondary_capture_ds)
        assert result.p_actual_length is None
        assert result.r_actual_length is None


class TestHologicArchivePRPresent:
    def test_hologic_archive_pr_present(self, hologic_archive_ds):
        """Full archive (P + R data) produces PR_PRESENT."""
        result = inspect_private_tags(hologic_archive_ds)
        assert result.diagnosis == PRDiagnosis.PR_PRESENT

    def test_hologic_archive_reports_private_creator(self, hologic_archive_ds):
        """Private creator string must be reported in the result."""
        result = inspect_private_tags(hologic_archive_ds)
        assert result.private_creator is not None
        assert "HOLOGIC" in result.private_creator.upper()

    def test_hologic_archive_actual_lengths_are_positive(self, hologic_archive_ds):
        """Actual P and R lengths must be positive."""
        result = inspect_private_tags(hologic_archive_ds)
        assert result.p_actual_length is not None and result.p_actual_length > 0
        assert result.r_actual_length is not None and result.r_actual_length > 0

    def test_hologic_archive_reports_sop_instance_uid(self, hologic_archive_ds):
        """SOP Instance UID must be captured in the result."""
        result = inspect_private_tags(hologic_archive_ds)
        assert result.sop_instance_uid not in ("", "UNKNOWN")


class TestHologicArchivePRPartial:
    def test_hologic_archive_pr_partial(self, hologic_archive_p_only_ds):
        """Archive with only P data (no R tags) produces PR_PARTIAL."""
        result = inspect_private_tags(hologic_archive_p_only_ds)
        assert result.diagnosis == PRDiagnosis.PR_PARTIAL

    def test_pr_partial_generates_warning(self, hologic_archive_p_only_ds):
        """PR_PARTIAL diagnosis must include at least one descriptive warning."""
        result = inspect_private_tags(hologic_archive_p_only_ds)
        assert len(result.warnings) > 0


class TestLengthMismatch:
    def test_length_mismatch_warning(self, hologic_archive_truncated_ds):
        """Declared length differs from actual by >1: length mismatch warning expected."""
        result = inspect_private_tags(hologic_archive_truncated_ds)
        mismatch_warnings = [w for w in result.warnings if "mismatch" in w.lower()]
        assert len(mismatch_warnings) > 0, (
            f"Expected a length-mismatch warning. Warnings were: {result.warnings}"
        )

    def test_length_mismatch_does_not_prevent_pr_present(
        self, hologic_archive_truncated_ds
    ):
        """A length mismatch is a warning, not a disqualification from PR_PRESENT."""
        result = inspect_private_tags(hologic_archive_truncated_ds)
        assert result.diagnosis == PRDiagnosis.PR_PRESENT


class TestPaddingByteAllowed:
    def test_padding_byte_allowed(self, hologic_archive_padded_ds):
        """Declared length = actual - 1 (DICOM odd-length padding): no length warning."""
        result = inspect_private_tags(hologic_archive_padded_ds)
        mismatch_warnings = [w for w in result.warnings if "mismatch" in w.lower()]
        assert len(mismatch_warnings) == 0, (
            f"Unexpected mismatch warning for padding byte. "
            f"Warnings: {result.warnings}"
        )

    def test_padding_byte_archive_is_pr_present(self, hologic_archive_padded_ds):
        """Padded archive with both blobs produces PR_PRESENT."""
        result = inspect_private_tags(hologic_archive_padded_ds)
        assert result.diagnosis == PRDiagnosis.PR_PRESENT


class TestEmptyBlobs:
    def test_empty_blobs_detected(self, hologic_archive_empty_ds):
        """Zero-length P/R blobs must NOT produce PR_PRESENT."""
        result = inspect_private_tags(hologic_archive_empty_ds)
        assert result.diagnosis != PRDiagnosis.PR_PRESENT, (
            "Empty blobs must not be treated as present P/R data."
        )

    def test_empty_blobs_actual_lengths_are_zero(self, hologic_archive_empty_ds):
        """Actual lengths for empty blobs should be zero."""
        result = inspect_private_tags(hologic_archive_empty_ds)
        assert result.p_actual_length == 0
        assert result.r_actual_length == 0


class TestNoPHI:
    def test_no_phi_in_result_from_secondary_capture(self, secondary_capture_ds):
        """Result from Secondary Capture must contain no patient-identifying fields."""
        result = inspect_private_tags(secondary_capture_ds)
        assert not _result_contains_phi(result)

    def test_no_phi_in_result_from_archive(self, hologic_archive_ds):
        """Result from Hologic archive must contain no patient-identifying fields."""
        result = inspect_private_tags(hologic_archive_ds)
        assert not _result_contains_phi(result)

    def test_result_does_not_store_binary_blobs(self, hologic_archive_ds):
        """raw_private_tags must store metadata only (length/VR), not raw bytes."""
        result = inspect_private_tags(hologic_archive_ds)
        for key, meta in result.raw_private_tags.items():
            assert "VR" in meta, f"Missing VR metadata for tag {key}"
            assert not isinstance(meta.get("value"), (bytes, bytearray)), (
                f"Raw bytes found in raw_private_tags[{key}]. "
                "Only length/VR metadata should be stored."
            )