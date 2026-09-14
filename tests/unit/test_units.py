"""Unit tests for hologic_dxa.maps.units."""
from __future__ import annotations

import numpy as np
import pytest

from hologic_dxa.maps.units import (
    KNOWN_AREAL_DENSITY_UNITS,
    convert_to_g_cm2,
    validate_areal_density_unit,
)


class TestValidateArealDensityUnit:
    @pytest.mark.parametrize("unit", ["g/cm2", "g/cm²", "mg/cm2", "mg/cm²"])
    def test_valid_units_pass(self, unit):
        validate_areal_density_unit(unit)  # must not raise

    @pytest.mark.parametrize("unit", ["kg/m2", "HU", "g/mm2", "", "g cm2"])
    def test_invalid_units_raise(self, unit):
        with pytest.raises(ValueError, match="Unrecognised"):
            validate_areal_density_unit(unit)

    def test_known_units_set_non_empty(self):
        assert len(KNOWN_AREAL_DENSITY_UNITS) >= 4


class TestConvertToGCm2:
    def test_g_cm2_identity(self):
        result = convert_to_g_cm2(1.5, "g/cm2")
        assert result == pytest.approx(1.5)

    def test_g_cm2_unicode_identity(self):
        result = convert_to_g_cm2(1.5, "g/cm²")
        assert result == pytest.approx(1.5)

    def test_mg_cm2_to_g_cm2(self):
        result = convert_to_g_cm2(1000.0, "mg/cm2")
        assert result == pytest.approx(1.0)

    def test_mg_cm2_unicode_to_g_cm2(self):
        result = convert_to_g_cm2(500.0, "mg/cm²")
        assert result == pytest.approx(0.5)

    def test_array_conversion(self):
        arr = np.array([1000.0, 2000.0, 500.0])
        result = convert_to_g_cm2(arr, "mg/cm2")
        expected = np.array([1.0, 2.0, 0.5])
        np.testing.assert_allclose(result, expected)

    def test_unknown_unit_raises(self):
        with pytest.raises(ValueError, match="Unrecognised"):
            convert_to_g_cm2(1.0, "kg/m2")

    def test_zero_value(self):
        result = convert_to_g_cm2(0.0, "g/cm2")
        assert result == pytest.approx(0.0)

    def test_small_value_precision(self):
        result = convert_to_g_cm2(1.0, "mg/cm2")
        assert result == pytest.approx(0.001, rel=1e-9)
