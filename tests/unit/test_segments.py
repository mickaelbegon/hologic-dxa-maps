"""Unit tests for hologic_dxa.segments (CoR loading, R-file parsing, slicing)."""
from __future__ import annotations

import json
import struct

import numpy as np
import pytest
from pydicom.dataset import Dataset

from hologic_dxa.segments import (
    APEX_SCALED_STATUS,
    CALIBRATION_STATUS,
    SEGMENT_DEFINITIONS,
    XIPHOID_FRACTION,
    SegmentProfile,
    SliceGeometry,
    TissueModel,
    bsp_segment_name,
    compute_segment_profiles,
    lateral_constraint,
    load_cor_json,
    load_whole_body_mu,
    parse_r_file_images,
    profiles_to_dataframe,
    rescale_to_apex,
    slice_segment,
    suggest_markers,
    to_calibration_bundle,
)
from hologic_dxa.segments.whole_body import PEDESTALS

ROW_M = SliceGeometry().row_size_m
DEFS = {d.name: d for d in SEGMENT_DEFINITIONS}


def _tlv(rtype: int, payload: bytes) -> bytes:
    return struct.pack("<HI", rtype, len(payload) + 6) + payload


def _r_file(planes: np.ndarray) -> bytes:
    n_planes, n_rows, _ = planes.shape
    descriptor = bytearray(26)
    struct.pack_into("<H", descriptor, 20, n_rows)
    struct.pack_into("<H", descriptor, 24, n_planes)
    body = planes.transpose(2, 1, 0).astype("<u2").tobytes()
    return _tlv(221, bytes(descriptor)) + _tlv(202, body) + _tlv(0, b"")


class TestLoadCorJson:
    def test_valid(self, tmp_path):
        p = tmp_path / "cor.json"
        p.write_text(json.dumps({"head_top": [50.0, 5.0], "wrist_L": [80, 90]}))
        cor = load_cor_json(p)
        assert cor["wrist_L"] == (80.0, 90.0)

    def test_unknown_marker(self, tmp_path):
        p = tmp_path / "cor.json"
        p.write_text(json.dumps({"nose": [1, 1]}))
        with pytest.raises(ValueError, match="unknown marker"):
            load_cor_json(p)

    def test_outside_grid(self, tmp_path):
        p = tmp_path / "cor.json"
        p.write_text(json.dumps({"head_top": [106.0, 5.0]}))
        with pytest.raises(ValueError, match="outside"):
            load_cor_json(p)

    def test_malformed_value(self, tmp_path):
        p = tmp_path / "cor.json"
        p.write_text(json.dumps({"head_top": [1.0]}))
        with pytest.raises(ValueError, match=r"\[row, col\]"):
            load_cor_json(p)


class TestRFile:
    def test_roundtrip_shape_and_values(self):
        rng = np.random.default_rng(0)
        planes = rng.integers(0, 4000, size=(6, 106, 150), dtype=np.uint16)
        images = parse_r_file_images(_r_file(planes))
        assert len(images) == 1
        np.testing.assert_array_equal(images[0], planes)

    def test_truncated_stream_returns_empty(self):
        assert parse_r_file_images(b"\x00\x01") == []

    def test_missing_tag_raises(self):
        with pytest.raises(ValueError, match="0023,1004"):
            load_whole_body_mu(Dataset())

    def test_pedestal_subtraction_and_energy_split(self):
        planes = np.zeros((6, 106, 150), dtype=np.uint16)
        for p, ped in enumerate(PEDESTALS):
            planes[p] = int(ped) + (100 if p % 2 == 0 else 40)
        ds = Dataset()
        ds.add_new((0x0023, 0x1004), "OB", _r_file(planes))
        mu = load_whole_body_mu(ds)
        assert mu.mu_h.shape == (106, 150)
        np.testing.assert_allclose(mu.mu_h, 100.0)
        np.testing.assert_allclose(mu.mu_l, 40.0)


def _block(rows: slice, cols: slice, value: float = 300.0) -> np.ndarray:
    img = np.zeros((106, 150))
    img[rows, cols] = value
    return img


class TestSlicing:
    def test_trunk_stops_at_shoulder_cor(self):
        mu_h = _block(slice(5, 100), slice(10, 90))
        cor = {
            "cervicothoracic": (53.0, 20.0),
            "lumbosacral": (53.0, 70.0),
            "glenohumeral_R": (30.0, 20.0),
            "glenohumeral_L": (76.0, 20.0),
        }
        geom = SliceGeometry()
        out = slice_segment(
            DEFS["trunk"], cor["cervicothoracic"], cor["lumbosacral"], mu_h, mu_h * 1.2,
            geometry=geom, allowed=lateral_constraint(DEFS["trunk"], cor, geom),
        )
        assert all(s is not None for s in out)
        for s in out:
            assert min(s.edge_rows) >= 30.0 - 0.5
            assert max(s.edge_rows) <= 76.0 + 0.5

    def test_thigh_does_not_cross_median_axis(self):
        mu_h = _block(slice(30, 77), slice(60, 140))
        cor = {
            "cervicothoracic": (53.0, 20.0),
            "lumbosacral": (53.0, 60.0),
            "hip_L": (60.0, 60.0),
            "knee_L": (60.0, 100.0),
            "hip_R": (46.0, 60.0),
            "knee_R": (46.0, 100.0),
        }
        geom = SliceGeometry()
        for name, sign in (("thigh_L", 1), ("thigh_R", -1)):
            d = DEFS[name]
            out = slice_segment(
                d, cor[d.proximal], cor[d.distal], mu_h, mu_h * 1.2,
                geometry=geom, allowed=lateral_constraint(d, cor, geom),
            )
            for s in out:
                assert s is not None
                assert all(sign * (r - 53.0) >= -0.5 for r in s.edge_rows)
                assert s.body_width_m < 0.15

    def test_without_constraint_thigh_spills(self):
        mu_h = _block(slice(30, 77), slice(60, 140))
        d = DEFS["thigh_L"]
        out = slice_segment(d, (60.0, 60.0), (60.0, 100.0), mu_h, mu_h * 1.2)
        assert min(min(s.edge_rows) for s in out if s) < 50.0

    def test_degenerate_axis_gives_none(self):
        mu_h = _block(slice(0, 106), slice(0, 150))
        out = slice_segment(DEFS["head"], (50.0, 50.0), (50.0, 50.5), mu_h, mu_h, n_slices=4)
        assert out == [None] * 4

    def test_forearm_needs_wrist(self):
        mu_h = _block(slice(0, 106), slice(0, 150))
        base = {"elbow_L": (80.0, 60.0)}
        names = [p.name for p in compute_segment_profiles(mu_h, mu_h, base)]
        assert "forearm_L" not in names
        with_wrist = {**base, "wrist_L": (85.0, 100.0)}
        names = [p.name for p in compute_segment_profiles(mu_h, mu_h, with_wrist)]
        assert names == ["forearm_L"]

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="same shape"):
            compute_segment_profiles(np.zeros((10, 10)), np.zeros((10, 11)), {})


class TestDataFrame:
    def test_flags_calibration_status(self):
        mu_h = _block(slice(40, 60), slice(10, 120))
        cor = {"hip_R": (50.0, 20.0), "knee_R": (50.0, 60.0)}
        df = profiles_to_dataframe(compute_segment_profiles(mu_h, mu_h * 1.2, cor, n_slices=5))
        assert len(df) == 5
        assert set(df["calibration_status"]) == {CALIBRATION_STATUS}
        assert df["density_kg_m3"].between(900, 1100).all()
        assert df["body_width_m"].iloc[0] == pytest.approx(20 * ROW_M, abs=0.01)


class TestNeckPelvisArmsHandsFeet:
    def test_neck_and_pelvis_segments_built(self):
        mu_h = _block(slice(0, 106), slice(0, 150))
        cor = {
            "atlanto_occipital": (53.0, 20.0),
            "cervicothoracic": (53.0, 30.0),
            "lumbosacral": (53.0, 80.0),
            "hip_L": (60.0, 100.0),
            "hip_R": (46.0, 100.0),
        }
        names = [p.name for p in compute_segment_profiles(mu_h, mu_h, cor)]
        assert names == ["neck", "trunk", "pelvis"]

    def test_pelvis_width_limited_around_hips(self):
        mu_h = _block(slice(0, 106), slice(0, 150))
        cor = {"lumbosacral": (53.0, 60.0), "hip_L": (60.0, 90.0), "hip_R": (46.0, 90.0)}
        geom = SliceGeometry()
        pelvis = compute_segment_profiles(mu_h, mu_h, cor, geometry=geom)[0]
        limit = 2 * (7 * ROW_M + geom.pelvis_margin_m)
        assert all(s.body_width_m <= limit + 0.01 for s in pelvis.slices if s)

    def test_upper_arm_medial_extent_limited_by_lateral_extent(self):
        # arm (rows 55-65) touches the trunk (rows 20-55): one contiguous body region
        mu_h = _block(slice(20, 66), slice(0, 100))
        cor = {"glenohumeral_L": (60.0, 20.0), "elbow_L": (60.0, 80.0)}
        d = DEFS["upper_arm_L"]
        with_rule = slice_segment(d, cor["glenohumeral_L"], cor["elbow_L"], mu_h, mu_h * 1.2)
        no_rule = slice_segment(
            d, cor["glenohumeral_L"], cor["elbow_L"], mu_h, mu_h * 1.2,
            geometry=SliceGeometry(symmetric_edge_tissues=frozenset()),
        )
        assert min(min(s.edge_rows) for s in with_rule if s) >= 54.0
        assert min(min(s.edge_rows) for s in no_rule if s) < 52.0

    def test_upper_arm_proximal_half_borrows_distal_density(self):
        mu_h = _block(slice(55, 66), slice(0, 100))
        cols = np.arange(150)[None, :]
        mu_l = mu_h * np.where(cols < 50, 2.0, 1.0)
        cor = {"glenohumeral_L": (60.0, 20.0), "elbow_L": (60.0, 80.0)}
        seg = compute_segment_profiles(mu_h, mu_l, cor, n_slices=10)[0]
        assert seg.name == "upper_arm_L"
        prox, dist = seg.slices[:5], seg.slices[5:]
        assert all(s.density_source == "borrowed_distal" for s in prox)
        assert all(s.density_source == "measured" for s in dist)
        assert {round(s.density_kg_m3, 6) for s in prox} == {round(dist[0].density_kg_m3, 6)}
        assert prox[0].r_mean == pytest.approx(dist[0].r_mean)

    def test_hand_and_foot_have_uniform_density(self):
        mu_h = _block(slice(55, 66), slice(0, 150))
        cols = np.arange(150)[None, :]
        mu_l = mu_h * (1.0 + cols / 100.0)
        cor = {
            "wrist_L": (60.0, 20.0), "hand_end_L": (60.0, 50.0),
            "ankle_L": (60.0, 90.0), "foot_end_L": (60.0, 120.0),
        }
        by_name = {p.name: p for p in compute_segment_profiles(mu_h, mu_l, cor)}
        assert set(by_name) == {"hand_L", "foot_L"}
        for seg in by_name.values():
            rho = {round(s.density_kg_m3, 6) for s in seg.slices if s}
            assert len(rho) == 1
            assert all(s.density_source == "uniform_segment_mean" for s in seg.slices if s)

    def test_new_markers_accepted_in_json(self, tmp_path):
        p = tmp_path / "cor.json"
        p.write_text(
            json.dumps(
                {"hand_end_L": [80, 90], "foot_end_R": [30, 140], "iliac_crest_L": [60, 80]}
            )
        )
        assert set(load_cor_json(p)) == {"hand_end_L", "foot_end_R", "iliac_crest_L"}


class TestSuggestMarkers:
    def test_xiphoid_on_the_c7_lumbosacral_line(self):
        cor = {"cervicothoracic": (50.0, 30.0), "lumbosacral": (56.0, 80.0)}
        row, col = suggest_markers(cor)["xiphoid"]
        assert col == pytest.approx(30.0 + XIPHOID_FRACTION * 50.0)
        assert row == pytest.approx(50.0 + XIPHOID_FRACTION * 6.0)

    def test_no_suggestion_without_spine_markers(self):
        assert suggest_markers({"cervicothoracic": (50.0, 30.0)}) == {}

    def test_placed_marker_is_not_suggested(self):
        cor = {
            "cervicothoracic": (50.0, 30.0),
            "lumbosacral": (56.0, 80.0),
            "xiphoid": (52.0, 55.0),
        }
        assert suggest_markers(cor) == {}


class TestBundleExport:
    def _limb_profiles(self):
        mu_h = _block(slice(40, 60), slice(0, 150))
        mu_l = mu_h * 1.1
        cor = {
            "hip_L": (50.0, 10.0), "knee_L": (50.0, 50.0), "ankle_L": (50.0, 90.0),
            "foot_end_L": (50.0, 120.0),
        }
        return compute_segment_profiles(mu_h, mu_l, cor, n_slices=6)

    def test_side_names_follow_bsp_convention(self):
        assert bsp_segment_name("thigh_L") == "thigh_left"
        assert bsp_segment_name("forearm_R") == "forearm_right"
        assert bsp_segment_name("trunk") == "trunk"

    def test_bundle_keys_and_regional_pooling(self):
        b = to_calibration_bundle(
            self._limb_profiles(), sex="male", provenance="unit", n_slices=6
        )
        assert b["calibration_status"] == CALIBRATION_STATUS
        assert b["slice_orientation"] == "proximal_to_distal"
        assert set(b["regional_fat_fraction"]) == {"thigh", "shank"}
        assert set(b["slice_fat_fraction"]) == {"thigh_left", "shank_left", "foot_left"}
        assert all(len(v) == 6 for v in b["slice_fat_fraction"].values())
        assert b["regional_fat_fraction"]["thigh"] == pytest.approx(
            b["slice_fat_fraction"]["thigh_left"][0]
        )

    def test_uniform_segment_repeats_its_value(self):
        b = to_calibration_bundle(
            self._limb_profiles(), sex="female", provenance="unit", n_slices=6
        )
        assert len(set(b["slice_fat_fraction"]["foot_left"])) == 1

    def test_segment_with_invalid_slice_is_omitted_not_padded(self):
        prof = self._limb_profiles()
        thigh = next(p for p in prof if p.name == "thigh_L")
        broken = SegmentProfile(
            thigh.name, thigh.tissue, thigh.proximal, thigh.distal, thigh.length_m,
            [None, *thigh.slices[1:]],
        )
        others = [p for p in prof if p.name != "thigh_L"]
        b = to_calibration_bundle([broken, *others], sex="male", provenance="u", n_slices=6)
        assert "thigh_left" not in b["slice_fat_fraction"]
        assert "1 of 6" in b["omitted_segments"]["thigh_left"]

    def test_invalid_sex_and_no_limb_raise(self):
        with pytest.raises(ValueError, match="sex"):
            to_calibration_bundle(self._limb_profiles(), sex="x", provenance="u", n_slices=6)
        with pytest.raises(ValueError, match="No limb"):
            to_calibration_bundle([], sex="male", provenance="u", n_slices=6)


class TestApexScaling:
    def _left_leg(self, ratio=1.1):
        mu_h = _block(slice(40, 60), slice(0, 150))
        cor = {
            "hip_L": (50.0, 10.0), "knee_L": (50.0, 50.0), "ankle_L": (50.0, 90.0),
            "foot_end_L": (50.0, 120.0),
        }
        return compute_segment_profiles(mu_h, mu_h * ratio, cor, n_slices=6)

    @staticmethod
    def _pooled(profiles):
        num = den = 0.0
        for seg in profiles:
            w = seg.length_m / len(seg.slices)
            for sl in seg.slices:
                if sl is not None:
                    num += sl.fat_fraction * sl.mu_h_integral * w
                    den += sl.mu_h_integral * w
        return num / den

    def test_region_reaches_apex_fraction(self):
        scaled, scalings, warns = rescale_to_apex(self._left_leg(), {"L Leg": (4000.0, 6000.0)})
        assert warns[:0] == []
        assert self._pooled(scaled) == pytest.approx(0.4, abs=1e-3)
        assert scalings[0].fat_fraction_apex == pytest.approx(0.4)
        assert scalings[0].factor == pytest.approx(scaled[0].slices[0].fat_scale)

    def test_shape_ratio_is_preserved_and_density_recomputed(self):
        original = self._left_leg(ratio=1.3)
        scaled, _, _ = rescale_to_apex(original, {"L Leg": (3000.0, 7000.0)})
        for a, b in zip(original, scaled, strict=True):
            fa = [s.fat_fraction for s in a.slices if s]
            fb = [s.fat_fraction for s in b.slices if s]
            assert all(x > 0 for x in fa)
            assert [y / x for x, y in zip(fa, fb, strict=True)] == pytest.approx(
                [fb[0] / fa[0]] * len(fa)
            )
        thigh = next(p for p in scaled if p.name == "thigh_L").slices[0]
        assert thigh.density_kg_m3 == pytest.approx(
            TissueModel().density_from_fat_fraction(thigh.fat_fraction, "thigh")
        )
        assert thigh.scaled_to_apex

    def test_inputs_are_not_modified(self):
        original = self._left_leg()
        before = original[0].slices[0].fat_fraction
        rescale_to_apex(original, {"L Leg": (4000.0, 6000.0)})
        assert original[0].slices[0].fat_fraction == before
        assert not original[0].slices[0].scaled_to_apex

    def test_incomplete_region_is_left_unscaled_with_warning(self):
        mu_h = _block(slice(40, 60), slice(0, 150))
        cor = {"hip_L": (50.0, 10.0), "knee_L": (50.0, 50.0)}
        profiles = compute_segment_profiles(mu_h, mu_h * 1.1, cor, n_slices=6)
        scaled, scalings, warns = rescale_to_apex(profiles, {"L Leg": (4000.0, 6000.0)})
        assert scalings == []
        assert any("segments missing" in w for w in warns)
        assert not scaled[0].slices[0].scaled_to_apex

    def test_missing_apex_values_warn(self):
        _, scalings, warns = rescale_to_apex(self._left_leg(), {})
        assert scalings == []
        assert any("no APEX" in w for w in warns)

    def test_unreachable_target_is_not_forced(self):
        # R = 2.2 gives a zero fat fraction everywhere: no factor can reach 40 %
        _, scalings, warns = rescale_to_apex(
            self._left_leg(ratio=2.2), {"L Leg": (4000.0, 6000.0)}
        )
        assert scalings == []
        assert any("unreachable" in w for w in warns)

    def test_high_target_saturates_without_exceeding_one(self):
        scaled, scalings, _ = rescale_to_apex(self._left_leg(), {"L Leg": (9000.0, 1000.0)})
        assert all(s.fat_fraction <= 1.0 for p in scaled for s in p.slices if s)
        assert self._pooled(scaled) == pytest.approx(0.9, abs=2e-3)
        assert scalings[0].factor > 1.0

    def test_bundle_status_reflects_scaling(self):
        scaled, _, _ = rescale_to_apex(self._left_leg(), {"L Leg": (4000.0, 6000.0)})
        bundle = to_calibration_bundle(scaled, sex="male", provenance="u", n_slices=6)
        assert bundle["calibration_status"] == APEX_SCALED_STATUS
        plain = to_calibration_bundle(self._left_leg(), sex="male", provenance="u", n_slices=6)
        assert plain["calibration_status"] == CALIBRATION_STATUS
