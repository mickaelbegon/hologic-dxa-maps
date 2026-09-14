"""Unit tests for hologic_dxa.maps.geometry."""
from __future__ import annotations

import numpy as np
import pytest

from hologic_dxa.maps.geometry import (
    pixel_area_cm2_from_spacing,
    pixel_coordinates_mm,
    validate_pixel_area,
)


class TestPixelAreaCm2FromSpacing:
    def test_square_pixel_1mm(self):
        area = pixel_area_cm2_from_spacing(1.0, 1.0)
        assert area == pytest.approx(0.01)  # (1mm/10)^2 = 0.01 cm²

    def test_rectangular_pixel(self):
        area = pixel_area_cm2_from_spacing(2.0, 5.0)
        assert area == pytest.approx(0.2 * 0.5)  # 0.1 cm²

    def test_zero_row_spacing_raises(self):
        with pytest.raises(ValueError, match="strictly positive"):
            pixel_area_cm2_from_spacing(0.0, 1.0)

    def test_zero_col_spacing_raises(self):
        with pytest.raises(ValueError, match="strictly positive"):
            pixel_area_cm2_from_spacing(1.0, 0.0)

    def test_negative_spacing_raises(self):
        with pytest.raises(ValueError):
            pixel_area_cm2_from_spacing(-1.0, 1.0)

    def test_large_spacing(self):
        area = pixel_area_cm2_from_spacing(10.0, 10.0)
        assert area == pytest.approx(1.0)  # 1 cm²


class TestPixelCoordinatesMm:
    def test_shape(self):
        x, y = pixel_coordinates_mm(4, 5, 1.0, 2.0)
        assert x.shape == (4, 5)
        assert y.shape == (4, 5)

    def test_origin_default(self):
        x, y = pixel_coordinates_mm(3, 4, 1.0, 1.0)
        assert x[0, 0] == pytest.approx(0.0)
        assert y[0, 0] == pytest.approx(0.0)

    def test_column_spacing_applied_to_x(self):
        x, y = pixel_coordinates_mm(2, 4, 1.0, 2.5)
        # x should increase by col_spacing along columns
        assert x[0, 1] - x[0, 0] == pytest.approx(2.5)
        assert x[0, 3] == pytest.approx(3 * 2.5)

    def test_row_spacing_applied_to_y(self):
        x, y = pixel_coordinates_mm(4, 2, 3.0, 1.0)
        # y should increase by row_spacing along rows
        assert y[1, 0] - y[0, 0] == pytest.approx(3.0)
        assert y[3, 0] == pytest.approx(3 * 3.0)

    def test_custom_origin(self):
        x, y = pixel_coordinates_mm(2, 2, 1.0, 1.0, origin_mm=(10.0, 20.0))
        assert x[0, 0] == pytest.approx(10.0)
        assert y[0, 0] == pytest.approx(20.0)
        assert x[0, 1] == pytest.approx(11.0)
        assert y[1, 0] == pytest.approx(21.0)

    def test_dtype_float64(self):
        x, y = pixel_coordinates_mm(2, 2, 1.0, 1.0)
        assert x.dtype == np.float64
        assert y.dtype == np.float64


class TestValidatePixelArea:
    def test_valid_scalar(self):
        validate_pixel_area(0.01)  # must not raise

    def test_valid_array(self):
        validate_pixel_area(np.full((3, 4), 0.05))  # must not raise

    def test_none_raises(self):
        with pytest.raises(ValueError, match="None"):
            validate_pixel_area(None)

    def test_zero_raises(self):
        with pytest.raises(ValueError, match="strictly positive"):
            validate_pixel_area(0.0)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="strictly positive"):
            validate_pixel_area(-0.01)

    def test_array_with_zero_raises(self):
        arr = np.array([0.01, 0.0, 0.02])
        with pytest.raises(ValueError, match="strictly positive"):
            validate_pixel_area(arr)

    def test_empty_array_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_pixel_area(np.array([]))
