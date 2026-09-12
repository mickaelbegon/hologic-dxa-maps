"""Unit tests for hologic_dxa.maps.integration.

This module is imported with ``pytest.importorskip`` so tests are
automatically skipped if the module does not exist yet.
"""

from __future__ import annotations

import numpy as np
import pytest

integration = pytest.importorskip(
    "hologic_dxa.maps.integration",
    reason="hologic_dxa.maps.integration is not yet implemented",
)

regional_mass_g = integration.regional_mass_g
center_of_mass_mm = integration.center_of_mass_mm
planar_moment_of_inertia = integration.planar_moment_of_inertia

try:
    from hypothesis import assume, given, settings
    from hypothesis import strategies as st
    import hypothesis.extra.numpy as hnp
    _HYPOTHESIS_AVAILABLE = True
except ImportError:
    _HYPOTHESIS_AVAILABLE = False


# ---------------------------------------------------------------------------
# regional_mass_g
# ---------------------------------------------------------------------------

class TestRegionalMassG:
    def test_regional_mass_g_simple(self):
        """3×3 all-ones array, pixel area 1 cm², full mask → mass = 9 g."""
        density = np.ones((3, 3), dtype=np.float64)
        mass = regional_mass_g(density, pixel_area_cm2=1.0)
        assert np.isclose(mass, 9.0), f"Expected 9.0 g, got {mass}"

    def test_regional_mass_g_with_mask(self):
        """Partial mask → mass equals sum over masked pixels only."""
        density = np.ones((3, 3), dtype=np.float64) * 2.0
        mask = np.zeros((3, 3), dtype=bool)
        mask[0, :] = True   # first row only → 3 pixels
        mass = regional_mass_g(density, pixel_area_cm2=1.0, mask=mask)
        expected = 3 * 2.0 * 1.0
        assert np.isclose(mass, expected), f"Expected {expected} g, got {mass}"

    def test_regional_mass_g_full_mask_matches_no_mask(self):
        """Explicitly all-True mask must give same result as no mask."""
        density = np.arange(9, dtype=np.float64).reshape(3, 3)
        mask = np.ones((3, 3), dtype=bool)
        mass_no_mask = regional_mass_g(density, pixel_area_cm2=0.5)
        mass_full_mask = regional_mass_g(density, pixel_area_cm2=0.5, mask=mask)
        assert np.isclose(mass_no_mask, mass_full_mask)

    def test_regional_mass_g_zero_density(self):
        """All-zero density → mass = 0."""
        density = np.zeros((5, 5), dtype=np.float64)
        mass = regional_mass_g(density, pixel_area_cm2=2.0)
        assert mass == 0.0

    def test_regional_mass_g_scales_with_pixel_area(self):
        """Doubling pixel area must double the mass."""
        density = np.ones((4, 4), dtype=np.float64)
        mass_1 = regional_mass_g(density, pixel_area_cm2=1.0)
        mass_2 = regional_mass_g(density, pixel_area_cm2=2.0)
        assert np.isclose(mass_2, 2.0 * mass_1)

    def test_regional_mass_g_raises_without_pixel_area(self):
        """Passing pixel_area_cm2=None must raise ValueError."""
        density = np.ones((3, 3), dtype=np.float64)
        with pytest.raises(ValueError, match="pixel_area"):
            regional_mass_g(density, pixel_area_cm2=None)

    def test_regional_mass_g_per_pixel_area_array(self):
        """Per-pixel area array must give element-wise product sum."""
        density = np.ones((2, 2), dtype=np.float64)
        area = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mass = regional_mass_g(density, pixel_area_cm2=area)
        expected = float((density * area).sum())
        assert np.isclose(mass, expected)


# ---------------------------------------------------------------------------
# center_of_mass_mm
# ---------------------------------------------------------------------------

class TestCenterOfMassMm:
    def test_center_of_mass_uniform_density(self):
        """Uniform areal density → CoM at geometric center of the grid."""
        rows, cols = 5, 7
        density = np.ones((rows, cols), dtype=np.float64)
        spacing = (1.0, 1.0)  # mm
        com = center_of_mass_mm(density, pixel_spacing_mm=spacing)
        # Geometric center for 0-indexed pixel positions: (rows-1)/2, (cols-1)/2
        expected_row = (rows - 1) / 2.0
        expected_col = (cols - 1) / 2.0
        assert len(com) == 2, "CoM must return a 2-tuple (row_mm, col_mm)"
        assert np.isclose(com[0], expected_row, atol=1e-9), (
            f"Row CoM: expected {expected_row}, got {com[0]}"
        )
        assert np.isclose(com[1], expected_col, atol=1e-9), (
            f"Col CoM: expected {expected_col}, got {com[1]}"
        )

    def test_center_of_mass_point_mass(self):
        """Single non-zero pixel → CoM at that pixel position in mm."""
        density = np.zeros((5, 5), dtype=np.float64)
        density[2, 3] = 10.0
        spacing = (2.0, 3.0)  # mm
        com = center_of_mass_mm(density, pixel_spacing_mm=spacing)
        expected_row_mm = 2 * spacing[0]
        expected_col_mm = 3 * spacing[1]
        assert np.isclose(com[0], expected_row_mm, atol=1e-9), (
            f"Row CoM: expected {expected_row_mm}, got {com[0]}"
        )
        assert np.isclose(com[1], expected_col_mm, atol=1e-9), (
            f"Col CoM: expected {expected_col_mm}, got {com[1]}"
        )

    def test_center_of_mass_with_mask(self):
        """Mask restricts the CoM computation to the masked region."""
        density = np.ones((5, 5), dtype=np.float64)
        mask = np.zeros((5, 5), dtype=bool)
        mask[4, 4] = True  # only bottom-right pixel
        com = center_of_mass_mm(density, pixel_spacing_mm=(1.0, 1.0), mask=mask)
        assert np.isclose(com[0], 4.0, atol=1e-9)
        assert np.isclose(com[1], 4.0, atol=1e-9)

    def test_center_of_mass_spacing_scales_result(self):
        """Doubling pixel spacing must double CoM coordinates."""
        density = np.ones((3, 3), dtype=np.float64)
        com1 = center_of_mass_mm(density, pixel_spacing_mm=(1.0, 1.0))
        com2 = center_of_mass_mm(density, pixel_spacing_mm=(2.0, 2.0))
        assert np.isclose(com2[0], 2.0 * com1[0], atol=1e-9)
        assert np.isclose(com2[1], 2.0 * com1[1], atol=1e-9)


# ---------------------------------------------------------------------------
# planar_moment_of_inertia
# ---------------------------------------------------------------------------

class TestPlanarMomentOfInertia:
    def test_moment_of_inertia_positive(self):
        """Any spatially distributed mass must have I_perp > 0."""
        density = np.ones((5, 5), dtype=np.float64)
        density[2, 2] = 0.0  # remove center to ensure non-trivial distribution
        I = planar_moment_of_inertia(
            density,
            pixel_area_cm2=1.0,
            pixel_spacing_mm=(1.0, 1.0),
        )
        assert I > 0.0, f"Expected positive moment of inertia, got {I}"

    def test_moment_of_inertia_point_at_origin(self):
        """Single pixel at the CoM (origin) → I_perp = 0."""
        density = np.zeros((5, 5), dtype=np.float64)
        density[0, 0] = 1.0  # single pixel at [0,0] → CoM at origin
        I = planar_moment_of_inertia(
            density,
            pixel_area_cm2=1.0,
            pixel_spacing_mm=(1.0, 1.0),
        )
        # CoM is at [0,0]; pixel is at [0,0]; distance = 0 → I = 0
        assert np.isclose(I, 0.0, atol=1e-12), (
            f"Single pixel at CoM: expected I=0, got {I}"
        )

    def test_moment_of_inertia_increases_with_spread(self):
        """Moving mass further from the centre must increase I_perp."""
        density_near = np.zeros((9, 9), dtype=np.float64)
        density_near[4, 4] = 1.0  # centre

        density_near[3, 4] = 1.0  # one pixel away
        density_near[5, 4] = 1.0

        density_far = np.zeros((9, 9), dtype=np.float64)
        density_far[4, 4] = 1.0  # centre

        density_far[0, 4] = 1.0  # four pixels away
        density_far[8, 4] = 1.0

        I_near = planar_moment_of_inertia(
            density_near, pixel_area_cm2=1.0, pixel_spacing_mm=(1.0, 1.0)
        )
        I_far = planar_moment_of_inertia(
            density_far, pixel_area_cm2=1.0, pixel_spacing_mm=(1.0, 1.0)
        )
        assert I_far > I_near, (
            f"Farther spread should give larger I: I_near={I_near}, I_far={I_far}"
        )

    def test_moment_of_inertia_non_negative(self):
        """Moment of inertia must never be negative."""
        rng = np.random.default_rng(seed=7)
        density = rng.uniform(0.0, 10.0, size=(8, 8))
        I = planar_moment_of_inertia(
            density,
            pixel_area_cm2=0.5,
            pixel_spacing_mm=(1.5, 1.5),
        )
        assert I >= 0.0


# ---------------------------------------------------------------------------
# Hypothesis property tests
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HYPOTHESIS_AVAILABLE, reason="hypothesis not installed")
class TestMassProperties:
    @given(
        st.arrays(
            np.float64,
            shape=(5, 5),
            elements=st.floats(0, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=100)
    def test_mass_invariant_to_array_transpose(self, areal_density, pixel_area):
        """Mass must be the same for an array and its transpose (same sum)."""
        mass_orig = regional_mass_g(areal_density, pixel_area_cm2=pixel_area)
        mass_transposed = regional_mass_g(
            areal_density.T, pixel_area_cm2=pixel_area
        )
        assert np.isclose(mass_orig, mass_transposed, rtol=1e-9, atol=0), (
            f"mass_orig={mass_orig}, mass_transposed={mass_transposed}"
        )

    @given(
        st.arrays(
            np.float64,
            shape=(5, 5),
            elements=st.floats(0, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=100)
    def test_mass_non_negative_for_non_negative_inputs(
        self, areal_density, pixel_area
    ):
        """Non-negative density × positive area → non-negative mass."""
        mass = regional_mass_g(areal_density, pixel_area_cm2=pixel_area)
        assert mass >= 0.0, f"Got negative mass {mass} for non-negative inputs"

    @given(
        st.arrays(
            np.float64,
            shape=(5, 5),
            elements=st.floats(0, 100, allow_nan=False, allow_infinity=False),
        ),
        st.floats(0.01, 10.0),
    )
    @settings(max_examples=50)
    def test_mass_scales_linearly_with_pixel_area(self, areal_density, pixel_area):
        """Doubling pixel area must exactly double mass (linearity)."""
        mass1 = regional_mass_g(areal_density, pixel_area_cm2=pixel_area)
        mass2 = regional_mass_g(areal_density, pixel_area_cm2=2.0 * pixel_area)
        assert np.isclose(mass2, 2.0 * mass1, rtol=1e-9, atol=1e-14), (
            f"mass1={mass1}, mass2={mass2}"
        )
