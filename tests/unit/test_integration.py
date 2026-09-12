"""Unit tests for hologic_dxa.maps.integration."""

from __future__ import annotations

import numpy as np
import pytest

integration = pytest.importorskip(
    "hologic_dxa.maps.integration",
    reason="hologic_dxa.maps.integration not available",
)
geometry_mod = pytest.importorskip("hologic_dxa.maps.geometry")

regional_mass_g = integration.regional_mass_g
center_of_mass_mm = integration.center_of_mass_mm
planar_moment_of_inertia = integration.planar_moment_of_inertia
pixel_coordinates_mm = geometry_mod.pixel_coordinates_mm

try:
    from hypothesis import given, settings
    from hypothesis import strategies as st
    import hypothesis.extra.numpy as hnp
    _HYPOTHESIS_AVAILABLE = True
except ImportError:
    _HYPOTHESIS_AVAILABLE = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _coords(rows: int, cols: int, spacing: tuple[float, float] = (1.0, 1.0)):
    return pixel_coordinates_mm(rows, cols, spacing[0], spacing[1])


def _full_mask(shape: tuple[int, int]) -> np.ndarray:
    return np.ones(shape, dtype=bool)


# ---------------------------------------------------------------------------
# regional_mass_g
# ---------------------------------------------------------------------------

class TestRegionalMassG:
    def test_simple_ones_full_mask(self):
        """3x3 all-ones, area=1 cm2, full mask -> mass = 9 g."""
        density = np.ones((3, 3), dtype=np.float64)
        mass = regional_mass_g(density, pixel_area_cm2=1.0, mask=_full_mask((3, 3)))
        assert np.isclose(mass, 9.0), f"Expected 9.0 g, got {mass}"

    def test_partial_mask(self):
        """Partial mask -> mass equals sum over masked pixels only."""
        density = np.ones((3, 3), dtype=np.float64) * 2.0
        mask = np.zeros((3, 3), dtype=bool)
        mask[0, :] = True  # first row: 3 pixels
        mass = regional_mass_g(density, pixel_area_cm2=1.0, mask=mask)
        assert np.isclose(mass, 6.0), f"Expected 6.0 g, got {mass}"

    def test_scales_with_pixel_area(self):
        """Doubling pixel area must double mass."""
        density = np.ones((4, 4), dtype=np.float64)
        mask = _full_mask((4, 4))
        m1 = regional_mass_g(density, pixel_area_cm2=1.0, mask=mask)
        m2 = regional_mass_g(density, pixel_area_cm2=2.0, mask=mask)
        assert np.isclose(m2, 2.0 * m1)

    def test_per_pixel_area_array(self):
        """Per-pixel area array gives element-wise product sum."""
        density = np.ones((2, 2), dtype=np.float64)
        area = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mask = _full_mask((2, 2))
        mass = regional_mass_g(density, pixel_area_cm2=area, mask=mask)
        assert np.isclose(mass, 10.0)

    def test_raises_without_pixel_area(self):
        """pixel_area_cm2=None must raise ValueError mentioning pixel_area."""
        density = np.ones((3, 3), dtype=np.float64)
        with pytest.raises((ValueError, TypeError)):
            regional_mass_g(density, pixel_area_cm2=None, mask=_full_mask((3, 3)))

    def test_empty_mask_raises(self):
        """All-False mask must raise ValueError."""
        density = np.ones((3, 3), dtype=np.float64)
        mask = np.zeros((3, 3), dtype=bool)
        with pytest.raises(ValueError):
            regional_mass_g(density, pixel_area_cm2=1.0, mask=mask)


# ---------------------------------------------------------------------------
# center_of_mass_mm
# ---------------------------------------------------------------------------

class TestCenterOfMassMm:
    def test_uniform_density_at_geometric_center(self):
        """Uniform density -> CoM at geometric center of the pixel grid."""
        rows, cols = 5, 7
        density = np.ones((rows, cols), dtype=np.float64)
        x_mm, y_mm = _coords(rows, cols, (1.0, 1.0))
        mask = _full_mask((rows, cols))
        cx, cy = center_of_mass_mm(density, 1.0, x_mm, y_mm, mask)
        assert np.isclose(cx, 3.0, atol=1e-9), f"cx={cx}"
        assert np.isclose(cy, 2.0, atol=1e-9), f"cy={cy}"

    def test_point_mass_at_known_position(self):
        """Single nonzero pixel -> CoM at that pixel's coordinates."""
        density = np.zeros((5, 5), dtype=np.float64)
        density[2, 3] = 10.0  # row=2, col=3
        x_mm, y_mm = _coords(5, 5, (2.0, 3.0))
        mask = np.zeros((5, 5), dtype=bool)
        mask[2, 3] = True
        cx, cy = center_of_mass_mm(density, 1.0, x_mm, y_mm, mask)
        assert np.isclose(cx, x_mm[2, 3], atol=1e-9)
        assert np.isclose(cy, y_mm[2, 3], atol=1e-9)

    def test_mask_restricts_com(self):
        """Mask on bottom-right pixel -> CoM there."""
        density = np.ones((5, 5), dtype=np.float64)
        x_mm, y_mm = _coords(5, 5, (1.0, 1.0))
        mask = np.zeros((5, 5), dtype=bool)
        mask[4, 4] = True
        cx, cy = center_of_mass_mm(density, 1.0, x_mm, y_mm, mask)
        assert np.isclose(cx, x_mm[4, 4], atol=1e-9)
        assert np.isclose(cy, y_mm[4, 4], atol=1e-9)


# ---------------------------------------------------------------------------
# planar_moment_of_inertia
# ---------------------------------------------------------------------------

class TestPlanarMomentOfInertia:
    def test_non_negative(self):
        """I_perp must always be >= 0."""
        rng = np.random.default_rng(42)
        density = rng.uniform(0.0, 10.0, (8, 8))
        x_mm, y_mm = _coords(8, 8, (1.5, 1.5))
        mask = _full_mask((8, 8))
        I = planar_moment_of_inertia(density, 0.5, x_mm, y_mm, mask)
        assert I >= 0.0, f"I_perp={I}"

    def test_point_at_com_gives_zero(self):
        """Single pixel at origin -> I_perp = 0 (distance = 0)."""
        density = np.zeros((5, 5), dtype=np.float64)
        density[0, 0] = 1.0
        x_mm, y_mm = _coords(5, 5, (1.0, 1.0))
        mask = np.zeros((5, 5), dtype=bool)
        mask[0, 0] = True
        I = planar_moment_of_inertia(density, 1.0, x_mm, y_mm, mask)
        assert np.isclose(I, 0.0, atol=1e-12), f"I={I}"

    def test_farther_spread_increases_inertia(self):
        """Mass farther from CoM -> larger I_perp."""
        x_mm, y_mm = _coords(9, 9, (1.0, 1.0))
        mask_all = _full_mask((9, 9))

        d_near = np.zeros((9, 9))
        d_near[4, 4] = 1.0
        d_near[3, 4] = 1.0
        d_near[5, 4] = 1.0

        d_far = np.zeros((9, 9))
        d_far[4, 4] = 1.0
        d_far[0, 4] = 1.0
        d_far[8, 4] = 1.0

        I_near = planar_moment_of_inertia(d_near, 1.0, x_mm, y_mm, mask_all)
        I_far = planar_moment_of_inertia(d_far, 1.0, x_mm, y_mm, mask_all)
        assert I_far > I_near, f"I_near={I_near}, I_far={I_far}"


# ---------------------------------------------------------------------------
# Hypothesis property tests
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HYPOTHESIS_AVAILABLE, reason="hypothesis not installed")
class TestMassProperties:
    @given(
        hnp.arrays(
            np.float64,
            shape=(5, 5),
            elements=st.floats(0, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=50)
    def test_mass_invariant_to_transpose(self, areal_density, pixel_area):
        """Mass must be the same for an array and its transpose."""
        mask = _full_mask((5, 5))
        mass_orig = regional_mass_g(areal_density, pixel_area, mask)
        mass_T = regional_mass_g(areal_density.T, pixel_area, _full_mask((5, 5)))
        assert np.isclose(mass_orig, mass_T, rtol=1e-9)

    @given(
        hnp.arrays(
            np.float64,
            shape=(5, 5),
            elements=st.floats(0, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=50)
    def test_mass_non_negative(self, areal_density, pixel_area):
        """Non-negative density x positive area -> non-negative mass."""
        mask = _full_mask((5, 5))
        mass = regional_mass_g(areal_density, pixel_area, mask)
        assert mass >= 0.0

    @given(
        hnp.arrays(
            np.float64,
            shape=st.tuples(st.integers(2, 5), st.integers(2, 5)),
            elements=st.floats(0.1, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=50)
    def test_com_within_coordinate_bounds(self, areal_density, pixel_area):
        """CoM must lie within the extent of the pixel coordinates."""
        rows, cols = areal_density.shape
        x_mm, y_mm = _coords(rows, cols, (1.0, 1.0))
        mask = _full_mask((rows, cols))
        cx, cy = center_of_mass_mm(areal_density, pixel_area, x_mm, y_mm, mask)
        assert x_mm.min() <= cx <= x_mm.max(), f"cx={cx} out of [{x_mm.min()}, {x_mm.max()}]"
        assert y_mm.min() <= cy <= y_mm.max(), f"cy={cy} out of [{y_mm.min()}, {y_mm.max()}]"

    @given(
        hnp.arrays(
            np.float64,
            shape=(4, 4),
            elements=st.floats(0.1, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=50)
    def test_inertia_non_negative(self, areal_density, pixel_area):
        """I_perp must never be negative."""
        x_mm, y_mm = _coords(4, 4, (1.0, 1.0))
        mask = _full_mask((4, 4))
        I = planar_moment_of_inertia(areal_density, pixel_area, x_mm, y_mm, mask)
        assert I >= 0.0, f"I_perp={I}"
