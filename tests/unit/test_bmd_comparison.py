"""Unit tests for hologic_dxa.dicom.bmd_comparison."""
from __future__ import annotations

import pytest
import pandas as pd

from hologic_dxa.dicom.bmd_comparison import (
    ComparisonResult,
    compare_apex_xml_vs_sr,
    comparison_summary,
    normalise_region,
    COMPARABLE_COLUMNS,
)
from hologic_dxa.dicom.hologic_xml import ApexRegionMeasurement, ApexResults, ApexScanInfo
from hologic_dxa.dicom.structured_report import SRMeasurement


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _xml(region: str, **kwargs) -> ApexRegionMeasurement:
    defaults = dict(area_cm2=None, bmc_g=None, bmd_g_cm2=None,
                    fat_g=None, lean_g=None, lean_bmc_g=None,
                    total_g=None, fat_percent=None, tscore=None, zscore=None)
    defaults.update(kwargs)
    return ApexRegionMeasurement(region=region, **defaults)


def _sr(code_value: str, code_meaning: str, value: float,
        region: str = "", unit_code: str = "g") -> SRMeasurement:
    return SRMeasurement(
        code_value=code_value,
        coding_scheme="DCM",
        code_meaning=code_meaning,
        value=value,
        unit_code=unit_code,
        unit_meaning=None,
        anatomic_region_code=None,
        anatomic_region_meaning=region,
        tracking_id=None,
        path_in_tree="root/0",
    )


def _make_xml_results(*measurements) -> ApexResults:
    return ApexResults(
        measurements=list(measurements),
        scan_info=ApexScanInfo(),
        warnings=[],
    )


# ---------------------------------------------------------------------------
# TestNormaliseRegion
# ---------------------------------------------------------------------------

class TestNormaliseRegion:
    @pytest.mark.parametrize("raw,expected", [
        ("L Arm",            "L Arm"),
        ("l arm",            "L Arm"),
        ("Left Arm",         "L Arm"),
        ("left arm",         "L Arm"),
        ("left upper extremity", "L Arm"),
        ("Right Arm",        "R Arm"),
        ("r arm",            "R Arm"),
        ("Lumbar Spine",     "L Spine"),
        ("lumbar",           "L Spine"),
        ("Thoracic Spine",   "T Spine"),
        ("Total Body",       "Total"),
        ("whole body",       "Total"),
        ("Subtotal",         "Subtotal"),
        ("sub total",        "Subtotal"),
        ("Unknown Region",   "Unknown Region"),   # passthrough
        ("  L Arm  ",        "L Arm"),            # whitespace stripped
    ])
    def test_normalise(self, raw, expected):
        assert normalise_region(raw) == expected


# ---------------------------------------------------------------------------
# TestCompareApexXmlVsSr — basic shape
# ---------------------------------------------------------------------------

class TestCompareApexXmlVsSrShape:
    def test_returns_comparison_result(self):
        xml = _make_xml_results(_xml("L Arm", bmc_g=147.1))
        sr  = [_sr("122490", "Bone Mineral Content", 147.0, region="Left Arm")]
        result = compare_apex_xml_vs_sr(xml, sr)
        assert isinstance(result, ComparisonResult)
        assert isinstance(result.df, pd.DataFrame)

    def test_df_has_required_columns(self):
        result = compare_apex_xml_vs_sr(_make_xml_results(), [])
        assert {"region", "quantity", "xml_value", "sr_value",
                "abs_diff", "rel_diff_pct"}.issubset(result.df.columns)

    def test_empty_sources_give_empty_df(self):
        result = compare_apex_xml_vs_sr(_make_xml_results(), [])
        assert len(result.df) == 0

    def test_xml_only_row_has_no_sr_value(self):
        xml = _make_xml_results(_xml("L Arm", bmc_g=147.1))
        result = compare_apex_xml_vs_sr(xml, [])
        row = result.df[result.df["quantity"] == "bmc_g"].iloc[0]
        assert row["xml_value"] == pytest.approx(147.1)
        assert pd.isna(row["sr_value"])

    def test_sr_only_row_has_no_xml_value(self):
        sr = [_sr("122490", "Bone Mineral Content", 147.0, region="Left Arm")]
        result = compare_apex_xml_vs_sr(_make_xml_results(), sr)
        row = result.df[result.df["quantity"] == "bmc_g"].iloc[0]
        assert pd.isna(row["xml_value"])
        assert row["sr_value"] == pytest.approx(147.0)


# ---------------------------------------------------------------------------
# TestCompareApexXmlVsSr — matched pairs
# ---------------------------------------------------------------------------

class TestCompareApexXmlVsSrMatched:
    @pytest.fixture
    def matched_result(self):
        xml = _make_xml_results(
            _xml("L Arm", bmc_g=147.1, bmd_g_cm2=0.954, area_cm2=154.3),
            _xml("R Arm", bmc_g=149.9, bmd_g_cm2=0.858),
        )
        sr = [
            _sr("122490", "Bone Mineral Content", 147.0, region="Left Arm"),
            _sr("122492", "Bone Mineral Density",   0.950, region="Left Arm"),
            _sr("122488", "Area",                  154.0, region="Left Arm"),
            _sr("122490", "Bone Mineral Content", 149.5, region="Right Arm"),
        ]
        return compare_apex_xml_vs_sr(xml, sr)

    def test_abs_diff_computed(self, matched_result):
        row = matched_result.df[
            (matched_result.df["region"] == "L Arm") &
            (matched_result.df["quantity"] == "bmc_g")
        ].iloc[0]
        assert row["abs_diff"] == pytest.approx(147.1 - 147.0)

    def test_rel_diff_pct_computed(self, matched_result):
        row = matched_result.df[
            (matched_result.df["region"] == "L Arm") &
            (matched_result.df["quantity"] == "bmc_g")
        ].iloc[0]
        expected_rel = 100.0 * (147.1 - 147.0) / 147.0
        assert row["rel_diff_pct"] == pytest.approx(expected_rel)

    def test_identical_values_give_zero_diff(self):
        xml = _make_xml_results(_xml("Total", bmc_g=2360.74))
        sr  = [_sr("122490", "Bone Mineral Content", 2360.74, region="Total")]
        result = compare_apex_xml_vs_sr(xml, sr)
        row = result.df[(result.df["region"] == "Total") &
                        (result.df["quantity"] == "bmc_g")].iloc[0]
        assert row["abs_diff"] == pytest.approx(0.0, abs=1e-9)
        assert row["rel_diff_pct"] == pytest.approx(0.0, abs=1e-9)

    def test_region_normalisation_matches_across_sources(self, matched_result):
        regions = matched_result.df["region"].unique()
        assert "L Arm" in regions
        assert "R Arm" in regions
        assert "Left Arm" not in regions   # should be normalised
        assert "Right Arm" not in regions

    def test_multiple_quantities_matched(self, matched_result):
        l_arm = matched_result.df[matched_result.df["region"] == "L Arm"]
        matched = l_arm.dropna(subset=["xml_value", "sr_value"])
        quantities = set(matched["quantity"])
        assert {"bmc_g", "bmd_g_cm2", "area_cm2"}.issubset(quantities)


# ---------------------------------------------------------------------------
# TestCompareApexXmlVsSr — unmatched reporting
# ---------------------------------------------------------------------------

class TestCompareApexXmlVsSrUnmatched:
    def test_unmatched_xml_regions_reported(self):
        xml = _make_xml_results(
            _xml("L Arm", bmc_g=147.0),
            _xml("Head",  bmc_g=590.0),
        )
        sr = [_sr("122490", "Bone Mineral Content", 147.0, region="Left Arm")]
        result = compare_apex_xml_vs_sr(xml, sr)
        assert "Head" in result.unmatched_xml_regions

    def test_unmatched_sr_regions_reported(self):
        xml = _make_xml_results(_xml("L Arm", bmc_g=147.0))
        sr = [
            _sr("122490", "Bone Mineral Content", 147.0, region="Left Arm"),
            _sr("122490", "Bone Mineral Content", 590.0, region="Head"),
        ]
        result = compare_apex_xml_vs_sr(xml, sr)
        assert "Head" in result.unmatched_sr_regions

    def test_unmapped_sr_codes_reported(self):
        xml = _make_xml_results(_xml("L Arm", bmc_g=147.0))
        sr = [_sr("999999", "Unknown Quantity", 1.0, region="Left Arm")]
        result = compare_apex_xml_vs_sr(xml, sr)
        codes = [c for c, _ in result.unmapped_sr_codes]
        assert "999999" in codes

    def test_sr_with_none_value_skipped(self):
        xml = _make_xml_results(_xml("L Arm", bmc_g=147.0))
        sr_m = _sr("122490", "Bone Mineral Content", 0.0, region="Left Arm")
        sr_m_none = SRMeasurement(
            code_value="122490", coding_scheme="DCM",
            code_meaning="Bone Mineral Content", value=None,
            unit_code="g", unit_meaning=None,
            anatomic_region_code=None, anatomic_region_meaning="Left Arm",
            tracking_id=None, path_in_tree="root/0",
        )
        result = compare_apex_xml_vs_sr(xml, [sr_m_none])
        # None-valued SR row should not appear as a matched pair
        row = result.df[(result.df["region"] == "L Arm") &
                        (result.df["quantity"] == "bmc_g")]
        assert len(row) == 1
        assert pd.isna(row.iloc[0]["sr_value"])


# ---------------------------------------------------------------------------
# TestComparisonSummary
# ---------------------------------------------------------------------------

class TestComparisonSummary:
    @pytest.fixture
    def result_with_matches(self):
        xml = _make_xml_results(
            _xml("L Arm", bmc_g=147.1, bmd_g_cm2=0.954),
            _xml("R Arm", bmc_g=149.9, bmd_g_cm2=0.858),
        )
        sr = [
            _sr("122490", "Bone Mineral Content", 147.0, region="Left Arm"),
            _sr("122492", "Bone Mineral Density",   0.950, region="Left Arm"),
            _sr("122490", "Bone Mineral Content", 150.0, region="Right Arm"),
            _sr("122492", "Bone Mineral Density",   0.860, region="Right Arm"),
        ]
        return compare_apex_xml_vs_sr(xml, sr)

    def test_returns_dataframe(self, result_with_matches):
        summary = comparison_summary(result_with_matches)
        assert isinstance(summary, pd.DataFrame)

    def test_has_required_columns(self, result_with_matches):
        summary = comparison_summary(result_with_matches)
        expected = {"quantity", "n_pairs", "mean_abs_diff", "max_abs_diff",
                    "mean_rel_diff_pct", "max_rel_diff_pct"}
        assert expected.issubset(summary.columns)

    def test_n_pairs_correct(self, result_with_matches):
        summary = comparison_summary(result_with_matches)
        bmc_row = summary[summary["quantity"] == "bmc_g"].iloc[0]
        assert bmc_row["n_pairs"] == 2

    def test_empty_result_gives_empty_summary(self):
        result = compare_apex_xml_vs_sr(_make_xml_results(), [])
        summary = comparison_summary(result)
        assert len(summary) == 0

    def test_mean_abs_diff_nonnegative(self, result_with_matches):
        summary = comparison_summary(result_with_matches)
        assert (summary["mean_abs_diff"] >= 0).all()

    def test_max_abs_diff_geq_mean(self, result_with_matches):
        summary = comparison_summary(result_with_matches)
        assert (summary["max_abs_diff"] >= summary["mean_abs_diff"]).all()
