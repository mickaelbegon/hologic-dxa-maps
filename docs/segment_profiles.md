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
