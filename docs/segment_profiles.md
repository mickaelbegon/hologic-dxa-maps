# Segment density profiles (experimental)

`hologic_dxa.segments` cuts a whole-body scan into body segments and returns a
density profile per segment. CLI: `hologic-dxa segment-profiles`.

**Status: experimental, not calibrated.** Fat fraction comes from an empirical
linear relation on `R = mu_L / mu_H`; bone fractions and tissue densities are
assumptions (the neck, pelvis, hand and foot values are placeholders to review).
Every output carries `calibration_status = empirical_uncalibrated`. Treat
results as relative profiles, not validated absolute densities.

## Input

- DICOM whole-body scan: `mu_H` and `mu_L` are rebuilt from the R file (tag
  `0023,1004`, see `dicom_hologic.md` section 8) after subtracting per-plane
  pedestals. Layout: axis 0 = lateral (row 0 = patient right), axis 1 =
  head-to-foot (col 0 = top of head), native pitch 5.7 mm x 13.1 mm.
- CoR JSON `{marker: [row, col]}` in that grid. Markers: `head_top`,
  `atlanto_occipital`, `cervicothoracic`, `lumbosacral`, and for each side
  (`_L`/`_R`, patient side) `glenohumeral`, `elbow`, `wrist`, `hand_end`, `hip`,
  `knee`, `ankle`, `foot_end`. Markers are placed manually; there is no automatic
  detection in the package. `hip_mid` (pelvis distal end) is derived from the hips.
  `xiphoid`, `iliac_crest_L/R` and `crotch` are accepted but not used yet (reserved
  for a thorax / abdomen split). The xiphoid is not visible on DXA: the editor
  suggests it at 48 % of the C7 to L5/S1 distance (about T9-T10, +-1 vertebra) for
  you to confirm or correct. The iliac crest replaces the umbilicus, which cannot
  be seen either (the crest lies 2-3 cm below its level).

## Segments

| Segment | Proximal | Distal |
|---|---|---|
| head | head_top | atlanto_occipital |
| neck | atlanto_occipital | cervicothoracic |
| trunk | cervicothoracic | lumbosacral |
| pelvis | lumbosacral | hip_mid |
| upper_arm, forearm | glenohumeral, elbow | elbow, wrist |
| hand | wrist | hand_end |
| thigh, shank | hip, knee | knee, ankle |
| foot | ankle | foot_end |

A segment is produced only when both bounding markers exist.

## Method

1. Each segment is cut into `n_slices` slices perpendicular to its CoR-to-CoR
   axis, computed in physical space because pixels are not square.
2. A slice keeps only the contiguous run of body pixels (`mu_H` > 50) that
   contains the axis, within a per-tissue half-width cap (`SliceGeometry`).
3. Lateral limits: trunk between the two glenohumeral CoR; pelvis within the hip
   half-span plus `pelvis_margin_m`; thigh on its own side of the median axis
   (line C7 to L5/S1, or the hip midpoint if the trunk CoR are missing).
4. Upper arm: the medial extent from the axis is limited to the lateral extent,
   so an arm held against the trunk cannot absorb trunk pixels.
5. Per slice: `R`, fat fraction, 2-component density `rho`, `mu_H` integral.
6. Upper arm: the proximal half of the slices takes its `R`/fat/density from the
   first distal-half slice (`density_source = borrowed_distal`); widths stay measured.
7. Hand and foot: one uniform density for the whole segment, from the pooled `R`
   weighted by the `mu_H` integral (`density_source = uniform_segment_mean`).

## Output

`segment_profiles.csv` (one row per valid slice) and `segment_profiles.json`
(parameters, provenance, SHA-256 of the input). File names never depend on DICOM
tag values.

## Known limits

- Forearm and shank have no medial/lateral rule other than the width cap; an
  adducted forearm touching the hip can still include hip pixels within the cap.
- The foot is seen in projection (toes point up when supine), so `foot_end` gives
  only the projected length.
- Pedestals and pixel pitch are empirical for Horizon W / APEX 13.6.

## BodyLoop calibration bundle

`hologic-dxa segment-profiles ... --sex male|female|other [--provenance LABEL]` also
writes `dxa_bundle.json`, read by `bodyloop_anthropometrics` (`--dxa`) as a
`DXACalibrationBundle`. The sex is never read from the DICOM: give it explicitly.

- `regional_fat_fraction` (Tier 1): limb types pooled over left and right, weighted by
  the `mu_H` integral of each slice.
- `slice_fat_fraction` (Tier 3): per BSP segment (`thigh_left`, ...), proximal to
  distal. A segment with any invalid slice is omitted and listed with its reason in
  `omitted_segments`; nothing is padded. Hand and foot repeat their uniform value.
- `segment_lengths_m`: CoR-to-CoR distances.
- `calibration_status = empirical_uncalibrated`: BodyLoop adds an `[A]` warning.
- `regional_bmc_g` is not produced, so Tier 2 cannot run from this bundle.

Compared with the Hologic APEX regional fat percentage (fat / total mass, arms and
legs, 6 scans), the slice-derived value is on average 2.9 points lower (SD 3.5); arms
differ more (mean -4.3) than legs (mean -1.4).

## Rescaling to Hologic APEX regional values

By default (`--apex-scale`, disable with `--no-apex-scale`) the limb fat fractions are
multiplied by one factor per region ("L Arm", "R Arm", "L Leg", "R Leg") so that their
mass-weighted mean equals the APEX `fat_g / (fat_g + lean_g)` (tag `0019,1000`). The
shape of the profile along the limb is kept; the level becomes Hologic's. The weight of a
slice is its `mu_H` integral times the slice thickness, a proxy for tissue mass.

- Assumes APEX "L Arm" is the patient's left arm and includes the hand (idem for legs
  and feet). Regions with a missing segment or missing APEX value are left unscaled,
  with a warning.
- Head, neck, trunk and pelvis are not rescaled.
- `calibration_status` becomes `apex_scaled_shape_uncalibrated` for rescaled rows and for
  a bundle whose limb segments are all rescaled; the profile shape stays empirical.
- Only fat is rescaled. APEX gives BMC per whole arm / leg, and splitting it among upper
  arm, forearm and hand would be a guess, so `regional_bmc_g` (Tier 2) is still absent.
- Factors on the 6 test scans: 0.96 to 1.31 (arms 0.96-1.31, legs 0.96-1.13); the
  uncalibrated relation underestimates arm fat most.
