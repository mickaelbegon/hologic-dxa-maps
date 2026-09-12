"""Unit tests for hologic_dxa.dicom.classify.classify."""

from __future__ import annotations

import pytest

from hologic_dxa.dicom.classify import DicomClass, QuantitativeEligibility, classify
from tests.synthetic_data.generators import (
    make_hologic_archive,
    make_parametric_map,
    make_secondary_capture,
    make_structured_report,
)


# ---------------------------------------------------------------------------
# Secondary Capture
# ---------------------------------------------------------------------------

class TestSecondaryCapture:
    def test_secondary_capture_not_eligible(self, secondary_capture_ds):
        """Secondary Capture → NOT_ELIGIBLE."""
        result = classify(secondary_capture_ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.NOT_ELIGIBLE

    def test_secondary_capture_class(self, secondary_capture_ds):
        """Secondary Capture → DICOM class is SECONDARY_CAPTURE."""
        result = classify(secondary_capture_ds)
        assert result.dicom_class == DicomClass.SECONDARY_CAPTURE

    def test_secondary_capture_cannot_be_used_as_quantitative_map(
        self, secondary_capture_ds
    ):
        """Diagnostic message must explicitly prohibit SC as a quantitative map.

        This is constraint #1 in docs/limitations.md: Secondary Capture pixel
        values are display-scaled and carry no physical unit. The diagnostic
        string must state this prohibition unambiguously.
        """
        result = classify(secondary_capture_ds)
        diagnostic_lower = result.diagnostic.lower()
        # The diagnostic must contain at least one of these prohibitive phrases
        prohibitive_keywords = ["must not", "cannot", "prohibited", "display-scaled"]
        has_prohibition = any(kw in diagnostic_lower for kw in prohibitive_keywords)
        assert has_prohibition, (
            f"Diagnostic message does not explicitly prohibit use of SC as a "
            f"quantitative map. Got: '{result.diagnostic}'"
        )


# ---------------------------------------------------------------------------
# Hologic Archive
# ---------------------------------------------------------------------------

class TestHologicArchive:
    def test_hologic_archive_classified_correctly(self, hologic_archive_ds):
        """Archive with P/R data → HOLOGIC_ARCHIVE class."""
        result = classify(hologic_archive_ds)
        assert result.dicom_class == DicomClass.HOLOGIC_ARCHIVE

    def test_hologic_archive_requires_calibration(self, hologic_archive_ds):
        """Archive → REQUIRES_CALIBRATION (cannot produce density without calibration)."""
        result = classify(hologic_archive_ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.REQUIRES_CALIBRATION

    def test_hologic_archive_diagnostic_mentions_calibration(self, hologic_archive_ds):
        """Diagnostic message for archive must mention calibration requirement."""
        result = classify(hologic_archive_ds)
        assert "calibration" in result.diagnostic.lower(), (
            f"Expected 'calibration' in diagnostic. Got: '{result.diagnostic}'"
        )


# ---------------------------------------------------------------------------
# Parametric Map
# ---------------------------------------------------------------------------

class TestParametricMap:
    def test_parametric_map_with_rwvm_eligible(self, parametric_map_ds):
        """Parametric Map with RWVM and float pixels → ELIGIBLE."""
        result = classify(parametric_map_ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.ELIGIBLE

    def test_parametric_map_dicom_class(self, parametric_map_ds):
        """Parametric Map → DICOM class is PARAMETRIC_MAP."""
        result = classify(parametric_map_ds)
        assert result.dicom_class == DicomClass.PARAMETRIC_MAP

    def test_parametric_map_has_float_pixel_flag(self, parametric_map_ds):
        """Parametric Map with FloatPixelData → has_float_pixel_data True."""
        result = classify(parametric_map_ds)
        assert result.has_float_pixel_data is True

    def test_parametric_map_has_rwvm_flag(self, parametric_map_ds):
        """Parametric Map with RWVM → has_real_world_value_mapping True."""
        result = classify(parametric_map_ds)
        assert result.has_real_world_value_mapping is True

    def test_parametric_map_without_rwvm_not_eligible(
        self, parametric_map_no_rwvm_ds
    ):
        """Parametric Map missing both RWVM and float pixels → NOT_ELIGIBLE."""
        result = classify(parametric_map_no_rwvm_ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.NOT_ELIGIBLE

    def test_parametric_map_without_rwvm_still_classified_as_pm(
        self, parametric_map_no_rwvm_ds
    ):
        """Parametric Map without RWVM is still classified as PARAMETRIC_MAP."""
        result = classify(parametric_map_no_rwvm_ds)
        assert result.dicom_class == DicomClass.PARAMETRIC_MAP

    def test_parametric_map_with_only_rwvm_eligible(self):
        """Parametric Map with RWVM but no float pixels → ELIGIBLE (RWVM is sufficient)."""
        ds = make_parametric_map(include_rwvm=True, include_float_pixels=False)
        result = classify(ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.ELIGIBLE

    def test_parametric_map_with_only_float_pixels_eligible(self):
        """Parametric Map with float pixels but no RWVM → ELIGIBLE (float pixels suffice)."""
        ds = make_parametric_map(include_rwvm=False, include_float_pixels=True)
        result = classify(ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.ELIGIBLE


# ---------------------------------------------------------------------------
# Structured Report
# ---------------------------------------------------------------------------

class TestStructuredReport:
    def test_sr_eligible(self, structured_report_ds):
        """SR → ELIGIBLE (may contain regional quantitative results)."""
        result = classify(structured_report_ds)
        assert result.quantitative_eligibility == QuantitativeEligibility.ELIGIBLE

    def test_sr_dicom_class(self, structured_report_ds):
        """SR → DICOM class is STRUCTURED_REPORT."""
        result = classify(structured_report_ds)
        assert result.dicom_class == DicomClass.STRUCTURED_REPORT


# ---------------------------------------------------------------------------
# Cross-cutting
# ---------------------------------------------------------------------------

class TestClassifyReturnsAllFields:
    def test_result_has_sop_class_uid(self, secondary_capture_ds):
        result = classify(secondary_capture_ds)
        assert result.sop_class_uid != ""

    def test_result_has_modality(self, secondary_capture_ds):
        result = classify(secondary_capture_ds)
        assert result.modality != ""

    def test_result_has_non_empty_diagnostic(self, secondary_capture_ds):
        result = classify(secondary_capture_ds)
        assert len(result.diagnostic) > 0
