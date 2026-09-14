"""Unit tests for hologic_dxa.dicom.structured_report."""
from __future__ import annotations

import pytest

from hologic_dxa.dicom.structured_report import (
    SRMeasurement,
    parse_structured_report,
    sr_to_dataframe,
)
from synthetic_data.generators import make_secondary_capture, make_structured_report


class TestParseStructuredReport:
    def test_returns_measurements_for_sr(self, structured_report_ds):
        result = parse_structured_report(structured_report_ds)
        assert isinstance(result, list)
        assert len(result) >= 1
        assert isinstance(result[0], SRMeasurement)

    def test_measurement_fields_populated(self, structured_report_ds):
        result = parse_structured_report(structured_report_ds)
        m = result[0]
        assert m.code_value == "122492"
        assert m.coding_scheme == "DCM"
        assert m.code_meaning == "Bone Mineral Density"
        assert m.value == pytest.approx(1.23)
        assert m.unit_code == "mg/cm2"

    def test_returns_empty_list_for_non_sr(self):
        ds = make_secondary_capture()
        result = parse_structured_report(ds)
        assert result == []

    def test_multiple_measurements(self):
        measurements_spec = [
            {
                "code_value": "122492",
                "coding_scheme": "DCM",
                "code_meaning": "Bone Mineral Density",
                "value": 1.23,
                "unit_code": "mg/cm2",
                "unit_meaning": "Milligrams per square centimeter",
            },
            {
                "code_value": "122489",
                "coding_scheme": "DCM",
                "code_meaning": "Bone Mineral Content",
                "value": 45.6,
                "unit_code": "g",
                "unit_meaning": "Gram",
            },
        ]
        ds = make_structured_report(measurements=measurements_spec)
        result = parse_structured_report(ds)
        assert len(result) == 2
        codes = {m.code_value for m in result}
        assert {"122492", "122489"} == codes

    def test_path_in_tree_populated(self, structured_report_ds):
        result = parse_structured_report(structured_report_ds)
        assert result[0].path_in_tree.startswith("root")

    def test_no_content_sequence_returns_empty(self):
        ds = make_structured_report()
        del ds.ContentSequence
        result = parse_structured_report(ds)
        assert result == []


class TestSrToDataframe:
    def test_columns_present(self, structured_report_ds):
        measurements = parse_structured_report(structured_report_ds)
        df = sr_to_dataframe(measurements)
        expected_cols = {
            "code_value", "coding_scheme", "code_meaning",
            "value", "unit_code", "unit_meaning",
            "anatomic_region_code", "anatomic_region_meaning",
            "tracking_id", "path_in_tree",
        }
        assert expected_cols.issubset(set(df.columns))

    def test_empty_list_returns_empty_dataframe(self):
        df = sr_to_dataframe([])
        assert len(df) == 0
        assert "code_value" in df.columns

    def test_row_count_matches_measurements(self, structured_report_ds):
        measurements = parse_structured_report(structured_report_ds)
        df = sr_to_dataframe(measurements)
        assert len(df) == len(measurements)

    def test_value_dtype_numeric(self, structured_report_ds):
        measurements = parse_structured_report(structured_report_ds)
        df = sr_to_dataframe(measurements)
        assert df["value"].dtype in (float, "float64", "object")
        assert float(df["value"].iloc[0]) == pytest.approx(1.23)
