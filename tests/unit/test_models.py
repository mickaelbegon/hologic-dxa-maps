"""Unit tests for hologic_dxa.models."""
from __future__ import annotations

import numpy as np
import pytest

from hologic_dxa.models import MapGeometry, QuantitativeMap, QuantitativeMapBundle
from hologic_dxa.provenance import Provenance
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_qmap(
    rows: int = 4,
    cols: int = 5,
    quantity: str = "fat_areal_density",
    unit: str = "g/cm2",
    pixel_area: float | None = 0.01,
) -> QuantitativeMap:
    rng = np.random.default_rng(0)
    values = rng.uniform(0.0, 1.0, size=(rows, cols))
    mask = np.ones((rows, cols), dtype=bool)
    return QuantitativeMap(
        values=values,
        quantity=quantity,
        unit=unit,
        valid_mask=mask,
        pixel_area_cm2=pixel_area,
        x_coordinates_mm=None,
        y_coordinates_mm=None,
        source_uid="1.2.3.4.5",
        provenance={"pipeline_version": "0.0.1"},
    )


def _make_geometry(rows: int = 4, cols: int = 5) -> MapGeometry:
    return MapGeometry(
        rows=rows,
        columns=cols,
        pixel_spacing_mm=(1.0, 1.0),
        origin_mm=(0.0, 0.0),
        is_fan_beam_corrected=False,
    )


def _make_provenance() -> Provenance:
    return Provenance(
        pipeline_version="0.0.1",
        timestamp_utc=datetime(2026, 1, 1, tzinfo=timezone.utc),
        input_hashes={},
        operations=[],
        source_sop_instance_uid="1.2.3",
        study_instance_uid="1.2.4",
        series_instance_uid="1.2.5",
        device_manufacturer="HOLOGIC, Inc.",
        device_model="Horizon W",
        apex_version="13.6",
        transfer_syntax_uid="1.2.840.10008.1.2.1",
    )


# ---------------------------------------------------------------------------
# QuantitativeMap
# ---------------------------------------------------------------------------

class TestQuantitativeMap:
    def test_construction_happy_path(self):
        qmap = _make_qmap()
        assert qmap.values.shape == (4, 5)
        assert qmap.values.dtype == np.float64

    def test_values_cast_to_float64(self):
        values = np.ones((3, 3), dtype=np.float32)
        mask = np.ones((3, 3), dtype=bool)
        qmap = QuantitativeMap(
            values=values,
            quantity="bmc",
            unit="g/cm2",
            valid_mask=mask,
            pixel_area_cm2=0.01,
            x_coordinates_mm=None,
            y_coordinates_mm=None,
            source_uid="1.2.3",
            provenance={},
        )
        assert qmap.values.dtype == np.float64

    def test_1d_values_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            QuantitativeMap(
                values=np.ones(10),
                quantity="bmc",
                unit="g/cm2",
                valid_mask=np.ones(10, dtype=bool),
                pixel_area_cm2=0.01,
                x_coordinates_mm=None,
                y_coordinates_mm=None,
                source_uid="1.2.3",
                provenance={},
            )

    def test_mask_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="valid_mask shape"):
            QuantitativeMap(
                values=np.ones((4, 5)),
                quantity="bmc",
                unit="g/cm2",
                valid_mask=np.ones((3, 5), dtype=bool),
                pixel_area_cm2=0.01,
                x_coordinates_mm=None,
                y_coordinates_mm=None,
                source_uid="1.2.3",
                provenance={},
            )

    def test_empty_quantity_raises(self):
        with pytest.raises(ValueError, match="quantity"):
            _make_qmap(quantity="")

    def test_empty_unit_raises(self):
        with pytest.raises(ValueError, match="unit"):
            _make_qmap(unit="")

    def test_empty_source_uid_raises(self):
        values = np.ones((2, 2))
        mask = np.ones((2, 2), dtype=bool)
        with pytest.raises(ValueError, match="source_uid"):
            QuantitativeMap(
                values=values,
                quantity="bmc",
                unit="g/cm2",
                valid_mask=mask,
                pixel_area_cm2=0.01,
                x_coordinates_mm=None,
                y_coordinates_mm=None,
                source_uid="",
                provenance={},
            )

    def test_pixel_area_none_raises_on_mass(self):
        qmap = _make_qmap(pixel_area=None)
        with pytest.raises(ValueError, match="pixel_area_cm2 is None"):
            qmap.integrated_mass_g()

    def test_integrated_mass_g_scalar_area(self):
        values = np.full((2, 2), 1.0)
        mask = np.ones((2, 2), dtype=bool)
        qmap = QuantitativeMap(
            values=values,
            quantity="fat",
            unit="g/cm2",
            valid_mask=mask,
            pixel_area_cm2=0.25,
            x_coordinates_mm=None,
            y_coordinates_mm=None,
            source_uid="1.2.3",
            provenance={},
        )
        assert qmap.integrated_mass_g() == pytest.approx(4 * 1.0 * 0.25)

    def test_integrated_mass_g_array_area(self):
        values = np.full((2, 2), 2.0)
        mask = np.ones((2, 2), dtype=bool)
        area = np.full((2, 2), 0.5)
        qmap = QuantitativeMap(
            values=values,
            quantity="lean",
            unit="g/cm2",
            valid_mask=mask,
            pixel_area_cm2=area,
            x_coordinates_mm=None,
            y_coordinates_mm=None,
            source_uid="1.2.3",
            provenance={},
        )
        assert qmap.integrated_mass_g() == pytest.approx(4 * 2.0 * 0.5)

    def test_pixel_count_valid(self):
        values = np.ones((3, 4))
        mask = np.array([
            [True, True, False, True],
            [False, True, True, True],
            [True, False, False, True],
        ])
        qmap = QuantitativeMap(
            values=values,
            quantity="fat",
            unit="g/cm2",
            valid_mask=mask,
            pixel_area_cm2=0.01,
            x_coordinates_mm=None,
            y_coordinates_mm=None,
            source_uid="1.2.3",
            provenance={},
        )
        assert qmap.pixel_count_valid() == int(mask.sum())


# ---------------------------------------------------------------------------
# MapGeometry
# ---------------------------------------------------------------------------

class TestMapGeometry:
    def test_construction(self):
        g = _make_geometry()
        assert g.rows == 4
        assert g.columns == 5

    def test_zero_rows_raises(self):
        with pytest.raises(ValueError, match="rows"):
            MapGeometry(
                rows=0, columns=5,
                pixel_spacing_mm=(1.0, 1.0), origin_mm=(0.0, 0.0),
                is_fan_beam_corrected=False,
            )

    def test_zero_columns_raises(self):
        with pytest.raises(ValueError, match="columns"):
            MapGeometry(
                rows=4, columns=0,
                pixel_spacing_mm=(1.0, 1.0), origin_mm=(0.0, 0.0),
                is_fan_beam_corrected=False,
            )

    def test_negative_spacing_raises(self):
        with pytest.raises(ValueError, match="positive"):
            MapGeometry(
                rows=4, columns=5,
                pixel_spacing_mm=(-1.0, 1.0), origin_mm=(0.0, 0.0),
                is_fan_beam_corrected=False,
            )

    def test_pixel_area_cm2_scalar(self):
        g = MapGeometry(
            rows=4, columns=5,
            pixel_spacing_mm=(2.0, 5.0), origin_mm=(0.0, 0.0),
            is_fan_beam_corrected=False,
        )
        assert g.pixel_area_cm2_scalar() == pytest.approx(0.2 * 0.5)

    def test_frozen_immutable(self):
        g = _make_geometry()
        with pytest.raises((TypeError, AttributeError)):
            g.rows = 99  # type: ignore[misc]


# ---------------------------------------------------------------------------
# QuantitativeMapBundle
# ---------------------------------------------------------------------------

class TestQuantitativeMapBundle:
    def _make_bundle(self, fat=True, lean=False, bmc=False, total=False) -> QuantitativeMapBundle:
        geometry = _make_geometry(4, 5)
        valid_mask = np.ones((4, 5), dtype=bool)
        fat_map = _make_qmap(4, 5) if fat else None
        lean_map = _make_qmap(4, 5, quantity="lean_areal_density") if lean else None
        bmc_map = _make_qmap(4, 5, quantity="bmc") if bmc else None
        total_map = _make_qmap(4, 5, quantity="total_areal_density") if total else None
        return QuantitativeMapBundle(
            fat=fat_map,
            lean=lean_map,
            bmc=bmc_map,
            total=total_map,
            body_mask=None,
            bone_mask=None,
            valid_mask=valid_mask,
            geometry=geometry,
            metadata={},
            provenance=_make_provenance(),
        )

    def test_at_least_one_map_required(self):
        geometry = _make_geometry()
        with pytest.raises(ValueError, match="at least one"):
            QuantitativeMapBundle(
                fat=None, lean=None, bmc=None, total=None,
                body_mask=None, bone_mask=None,
                valid_mask=np.ones((4, 5), dtype=bool),
                geometry=geometry, metadata={},
                provenance=_make_provenance(),
            )

    def test_construction_fat_only(self):
        bundle = self._make_bundle(fat=True)
        assert bundle.fat is not None
        assert bundle.lean is None

    def test_shape_mismatch_raises(self):
        geometry = _make_geometry(4, 5)
        wrong_shape_map = _make_qmap(3, 5)
        with pytest.raises(ValueError, match="shape"):
            QuantitativeMapBundle(
                fat=wrong_shape_map, lean=None, bmc=None, total=None,
                body_mask=None, bone_mask=None,
                valid_mask=np.ones((4, 5), dtype=bool),
                geometry=geometry, metadata={},
                provenance=_make_provenance(),
            )

    def test_validate_bundle_no_warnings_all_positive(self):
        bundle = self._make_bundle(fat=True)
        warnings = bundle.validate_bundle()
        # Warnings about missing pixel_area are acceptable; no shape/negative errors
        shape_warnings = [w for w in warnings if "Shape" in w or "negative" in w]
        assert shape_warnings == []

    def test_validate_bundle_warns_no_pixel_area(self):
        geometry = _make_geometry()
        fat_map = _make_qmap(4, 5, pixel_area=None)
        bundle = QuantitativeMapBundle(
            fat=fat_map, lean=None, bmc=None, total=None,
            body_mask=None, bone_mask=None,
            valid_mask=np.ones((4, 5), dtype=bool),
            geometry=geometry, metadata={},
            provenance=_make_provenance(),
        )
        warnings = bundle.validate_bundle()
        assert any("pixel_area" in w for w in warnings)

    def test_validate_bundle_warns_negative_values(self):
        geometry = _make_geometry(2, 2)
        values = np.array([[-0.5, 0.1], [0.2, 0.3]])
        mask = np.ones((2, 2), dtype=bool)
        fat_map = QuantitativeMap(
            values=values, quantity="fat", unit="g/cm2",
            valid_mask=mask, pixel_area_cm2=0.01,
            x_coordinates_mm=None, y_coordinates_mm=None,
            source_uid="1.2.3", provenance={},
        )
        bundle = QuantitativeMapBundle(
            fat=fat_map, lean=None, bmc=None, total=None,
            body_mask=None, bone_mask=None,
            valid_mask=mask, geometry=geometry, metadata={},
            provenance=_make_provenance(),
        )
        warnings = bundle.validate_bundle()
        assert any("negative" in w.lower() for w in warnings)

    def test_validate_bundle_warns_total_mismatch(self):
        geometry = _make_geometry(2, 2)
        mask = np.ones((2, 2), dtype=bool)
        fat_vals = np.full((2, 2), 0.5)
        lean_vals = np.full((2, 2), 0.3)
        total_vals = np.full((2, 2), 1.5)  # != fat + lean (0.8)

        def _qm(vals, qty):
            return QuantitativeMap(
                values=vals, quantity=qty, unit="g/cm2",
                valid_mask=mask, pixel_area_cm2=0.01,
                x_coordinates_mm=None, y_coordinates_mm=None,
                source_uid="1.2.3", provenance={},
            )

        bundle = QuantitativeMapBundle(
            fat=_qm(fat_vals, "fat"),
            lean=_qm(lean_vals, "lean"),
            bmc=None,
            total=_qm(total_vals, "total"),
            body_mask=None, bone_mask=None,
            valid_mask=mask, geometry=geometry, metadata={},
            provenance=_make_provenance(),
        )
        warnings = bundle.validate_bundle()
        assert any("Total" in w or "fat+lean" in w for w in warnings)
