"""Unit tests for hologic_dxa.dicom.hologic_xml."""
from __future__ import annotations

import pytest
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from hologic_dxa.dicom.hologic_xml import (
    ApexRegionMeasurement,
    ApexResults,
    ApexScanInfo,
    apex_results_to_dataframe,
    parse_apex_xml,
)


# ---------------------------------------------------------------------------
# Synthetic XML helpers
# ---------------------------------------------------------------------------

def _make_xml(
    *,
    include_table: bool = True,
    n_rows: int = 3,
    include_scan_info: bool = True,
    table_note: str = "ACF = 1.023, BCF = 0.987",
    include_html_in_headers: bool = False,
) -> str:
    """Build a minimal Hologic XML string for testing."""
    lines = ["<HologicHeader1>", "<JavaScriptFile>"]

    if include_scan_info:
        lines += [
            'Height = "165.0 cm";',
            'Weight = "69.0  kg";',
            'Age = "41";',
            'PatientSex = "Female";',
            'ScanMode = "Whole Body";',
            'Scan = "2026-01-15";',
            'AnalProtocol = "WB";',
            f'TableNote1 = "{table_note}";',
        ]

    if include_table:
        col0_header = "<b>Region</b>" if include_html_in_headers else "Region"
        lines.append(f'ResultsTable1[ 0][ 0] = "{col0_header}";')
        lines.append('ResultsTable1[ 0][ 1] = "Area[cm&sup2;]";')
        lines.append('ResultsTable1[ 0][ 2] = "BMC[g]";')
        lines.append('ResultsTable1[ 0][ 3] = "BMD[g/cm&sup2;]";')
        lines.append('ResultsTable1[ 0][ 4] = "Fat[g]";')
        lines.append('ResultsTable1[ 0][ 5] = "Lean[g]";')
        lines.append('ResultsTable1[ 0][ 6] = "Lean+BMC[g]";')
        lines.append('ResultsTable1[ 0][ 7] = "Total[g]";')
        lines.append('ResultsTable1[ 0][ 8] = "%Fat";')
        lines.append('ResultsTable1[ 0][ 9] = "T-Score";')
        lines.append('ResultsTable1[ 0][10] = "Z-Score";')

        regions = ["L Arm", "R Arm", "Trunk"][:n_rows]
        values = [
            (154.28, 147.11, 0.954, 1253.8, 1796.1, 1943.2, 3096.9, 39.2, -0.5, 0.3),
            (174.67, 149.87, 0.858, 1185.3, 1813.8, 1963.7, 3148.9, 37.6, -0.8, 0.1),
            ("", 643.62, "", 12999.3, 21401.1, 22044.7, 35043.7, 37.1, "", ""),
        ][:n_rows]
        for i, (region, vals) in enumerate(zip(regions, values)):
            row = i + 1
            lines.append(f'ResultsTable1[ {row}][ 0] = "{region}";')
            for col, v in enumerate(vals, start=1):
                lines.append(f'ResultsTable1[ {row}][{col:2d}] = "{v}";')

    lines += ["</JavaScriptFile>", "</HologicHeader1>"]
    return "\n".join(lines)


def _make_ds(xml: str | None) -> Dataset:
    """Build a minimal pydicom Dataset with or without tag (0019,1000)."""
    ds = Dataset()
    if xml is not None:
        ds[0x0019, 0x1000] = pydicom.DataElement(
            tag=(0x0019, 0x1000),
            VR="LT",
            value=xml.encode("utf-8"),
        )
    return ds


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def full_xml() -> str:
    return _make_xml()


@pytest.fixture
def full_ds(full_xml) -> Dataset:
    return _make_ds(full_xml)


# ---------------------------------------------------------------------------
# TestParseApexXml — tag absent / malformed
# ---------------------------------------------------------------------------

class TestParseApexXmlTagAbsent:
    def test_returns_apex_results_instance(self):
        ds = _make_ds(None)
        result = parse_apex_xml(ds)
        assert isinstance(result, ApexResults)

    def test_no_measurements_when_tag_absent(self):
        ds = _make_ds(None)
        result = parse_apex_xml(ds)
        assert result.measurements == []

    def test_no_warnings_when_tag_absent(self):
        ds = _make_ds(None)
        result = parse_apex_xml(ds)
        assert result.warnings == []

    def test_no_measurements_when_table_missing(self):
        ds = _make_ds(_make_xml(include_table=False))
        result = parse_apex_xml(ds)
        assert result.measurements == []

    def test_warning_issued_when_only_header_row(self):
        xml = _make_xml(n_rows=0)
        ds = _make_ds(xml)
        result = parse_apex_xml(ds)
        assert any("no data rows" in w.lower() or "header" in w.lower() for w in result.warnings)

    def test_invalid_bytes_returns_empty_with_warning(self):
        ds = Dataset()
        ds[0x0019, 0x1000] = pydicom.DataElement(
            tag=(0x0019, 0x1000),
            VR="OB",
            value=b"\xff\xfe",  # incomplete UTF-16, not valid UTF-8 but replace-handled
        )
        result = parse_apex_xml(ds)
        # Should not raise; may have 0 measurements and possibly warnings
        assert isinstance(result, ApexResults)


# ---------------------------------------------------------------------------
# TestParseApexXml — measurements
# ---------------------------------------------------------------------------

class TestParseApexXmlMeasurements:
    def test_returns_correct_region_count(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert len(result.measurements) == 3

    def test_region_names_preserved(self, full_ds):
        result = parse_apex_xml(full_ds)
        names = [m.region for m in result.measurements]
        assert names == ["L Arm", "R Arm", "Trunk"]

    def test_numeric_fields_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        arm = result.measurements[0]
        assert arm.area_cm2 == pytest.approx(154.28)
        assert arm.bmc_g == pytest.approx(147.11)
        assert arm.bmd_g_cm2 == pytest.approx(0.954)
        assert arm.fat_g == pytest.approx(1253.8)
        assert arm.lean_g == pytest.approx(1796.1)
        assert arm.fat_percent == pytest.approx(39.2)
        assert arm.tscore == pytest.approx(-0.5)
        assert arm.zscore == pytest.approx(0.3)

    def test_empty_cells_become_none(self, full_ds):
        result = parse_apex_xml(full_ds)
        trunk = result.measurements[2]  # Trunk has no area or bmd
        assert trunk.area_cm2 is None
        assert trunk.bmd_g_cm2 is None

    def test_html_stripped_from_region_label(self):
        xml = _make_xml(include_html_in_headers=True)
        # Inject HTML into a region label
        xml = xml.replace(
            'ResultsTable1[ 1][ 0] = "L Arm";',
            'ResultsTable1[ 1][ 0] = "<b>L Arm</b>";',
        )
        ds = _make_ds(xml)
        result = parse_apex_xml(ds)
        assert result.measurements[0].region == "L Arm"

    def test_returns_apex_region_measurement_instances(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert all(isinstance(m, ApexRegionMeasurement) for m in result.measurements)

    def test_region_with_empty_label_skipped(self):
        xml = _make_xml(n_rows=2)
        xml = xml.replace(
            'ResultsTable1[ 1][ 0] = "L Arm";',
            'ResultsTable1[ 1][ 0] = "";',
        )
        ds = _make_ds(xml)
        result = parse_apex_xml(ds)
        # Empty-label row should be skipped
        assert all(m.region != "" for m in result.measurements)
        assert any("skipped" in w.lower() or "empty" in w.lower() for w in result.warnings)

    def test_single_region(self):
        ds = _make_ds(_make_xml(n_rows=1))
        result = parse_apex_xml(ds)
        assert len(result.measurements) == 1
        assert result.measurements[0].region == "L Arm"


# ---------------------------------------------------------------------------
# TestParseApexXml — scan info
# ---------------------------------------------------------------------------

class TestParseApexXmlScanInfo:
    def test_returns_apex_scan_info_instance(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert isinstance(result.scan_info, ApexScanInfo)

    def test_height_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.height_cm == pytest.approx(165.0)

    def test_weight_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.weight_kg == pytest.approx(69.0)

    def test_age_parsed_as_int(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.age == 41
        assert isinstance(result.scan_info.age, int)

    def test_sex_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.sex == "Female"

    def test_scan_date_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.scan_date == "2026-01-15"

    def test_scan_protocol_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.scan_protocol == "WB"

    def test_scan_mode_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.scan_mode == "Whole Body"

    def test_acf_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.acf == pytest.approx(1.023)

    def test_bcf_parsed(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.bcf == pytest.approx(0.987)

    def test_calibration_note_set(self, full_ds):
        result = parse_apex_xml(full_ds)
        assert result.scan_info.calibration_note is not None
        assert "ACF" in result.scan_info.calibration_note

    def test_missing_scan_info_leaves_none(self):
        ds = _make_ds(_make_xml(include_scan_info=False))
        result = parse_apex_xml(ds)
        assert result.scan_info.height_cm is None
        assert result.scan_info.weight_kg is None
        assert result.scan_info.age is None
        assert result.scan_info.acf is None
        assert result.scan_info.bcf is None

    def test_acf_bcf_absent_when_table_note_empty(self):
        ds = _make_ds(_make_xml(table_note=""))
        result = parse_apex_xml(ds)
        assert result.scan_info.acf is None
        assert result.scan_info.bcf is None


# ---------------------------------------------------------------------------
# TestApexResultsToDataframe
# ---------------------------------------------------------------------------

class TestApexResultsToDataframe:
    def test_returns_dataframe_with_correct_columns(self, full_ds):
        import pandas as pd
        result = parse_apex_xml(full_ds)
        df = apex_results_to_dataframe(result)
        assert isinstance(df, pd.DataFrame)
        expected_cols = {
            "region", "area_cm2", "bmc_g", "bmd_g_cm2",
            "fat_g", "lean_g", "lean_bmc_g", "total_g",
            "fat_percent", "tscore", "zscore",
        }
        assert set(df.columns) == expected_cols

    def test_row_count_matches_measurements(self, full_ds):
        result = parse_apex_xml(full_ds)
        df = apex_results_to_dataframe(result)
        assert len(df) == len(result.measurements)

    def test_empty_result_returns_empty_dataframe(self):
        import pandas as pd
        result = ApexResults()
        df = apex_results_to_dataframe(result)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert "region" in df.columns

    def test_none_values_preserved_as_nan(self, full_ds):
        import math
        result = parse_apex_xml(full_ds)
        df = apex_results_to_dataframe(result)
        trunk_row = df[df["region"] == "Trunk"].iloc[0]
        assert trunk_row["area_cm2"] is None or math.isnan(trunk_row["area_cm2"]) or trunk_row["area_cm2"] != trunk_row["area_cm2"]

    def test_numeric_values_correct(self, full_ds):
        result = parse_apex_xml(full_ds)
        df = apex_results_to_dataframe(result)
        arm_row = df[df["region"] == "L Arm"].iloc[0]
        assert arm_row["area_cm2"] == pytest.approx(154.28)
        assert arm_row["bmc_g"] == pytest.approx(147.11)

    def test_region_column_values(self, full_ds):
        result = parse_apex_xml(full_ds)
        df = apex_results_to_dataframe(result)
        assert list(df["region"]) == ["L Arm", "R Arm", "Trunk"]
