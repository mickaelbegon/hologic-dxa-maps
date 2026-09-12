# Pipeline Limitations

This document describes the fundamental limitations of the `hologic-dxa-maps`
pipeline. Each entry is marked **confirmed**, **suspected**, or **unknown**
depending on the level of evidence available.

---

## 1. Secondary Capture Images Are Prohibited as Quantitative Maps

**Confirmed — non-negotiable constraint.**

Secondary Capture (SC) DICOM objects (SOP Class 1.2.840.10008.5.1.4.1.1.7)
are display-rendered images. Their pixel values are:

- Windowed and level-shifted for visual presentation.
- Clipped to the display bit depth (typically 8-bit after windowing).
- Irreversibly processed: the mapping from raw detector signal to display value
  is lossy and not invertible without the original calibration look-up tables.

**Secondary Capture images MUST NOT be used as quantitative maps under any
circumstances.** The pipeline enforces this constraint in
`hologic_dxa.dicom.classify`: any dataset whose DICOM class is
`SECONDARY_CAPTURE` receives `QuantitativeEligibility.NOT_ELIGIBLE`.
Overriding this classification is not supported.

---

## 2. Proprietary P/R File Format — Cannot Decode Without SDK

**Confirmed limitation.**

Hologic DXA archive objects embed two binary files:

- **P file** — believed to contain projection data or processed detector counts.
- **R file** — believed to contain reference or calibration-related data.

These files are stored in private DICOM tags (group `0023`) using a proprietary
format. Hologic has not published documentation for the internal structure of
these files. Without the official Hologic SDK or a reverse-engineering effort,
decoding the binary content is not possible.

The pipeline can detect and report the presence of P/R data (`PRDiagnosis.PR_PRESENT`)
but **cannot read or interpret the embedded binary content**. Any attempt to
treat P/R blobs as raw areal density values would be incorrect.

---

## 3. Absence of Calibration — Reconstruction Is Blocked

**Confirmed limitation.**

Even if the P/R file format were known, converting raw detector signal to
calibrated areal density (g/cm²) requires:

1. **Bone-equivalent phantom calibration factors** — acquired from periodic
   phantom scans, proprietary to each scanner installation.
2. **Beam-hardening correction tables** — specific to the X-ray tube spectrum
   and filtration in use.
3. **Detector linearity corrections** — factory-calibrated and stored in the
   scanner firmware.

None of these factors are available through standard DICOM tags. They reside in
the scanner's internal database and are not exported alongside patient DICOMs.

**Consequence:** the pipeline classifies Hologic archive objects as
`REQUIRES_CALIBRATION` and cannot proceed to quantitative reconstruction without
external supply of these factors. See `docs/dicom_hologic.md` for what
information can be read from DICOM metadata.

---

## 4. Fundamental Underdetermination: 2 Measurements, 3 Tissue Components

**Confirmed — fundamental physics constraint.**

Dual-energy DXA uses two X-ray energies (Low, L and High, H) to acquire two
attenuation measurements per pixel:

$$a_L = \sigma_{\text{fat}} \cdot \rho_{\text{fat}} + \sigma_{\text{lean}} \cdot \rho_{\text{lean}} + \sigma_{\text{BMC}} \cdot \rho_{\text{BMC}}$$

$$a_H = \sigma_{\text{fat}}' \cdot \rho_{\text{fat}} + \sigma_{\text{lean}}' \cdot \rho_{\text{lean}} + \sigma_{\text{BMC}}' \cdot \rho_{\text{BMC}}$$

where $\sigma$ are the known mass-attenuation coefficients and $\rho$ are the
unknown areal densities (g/cm²) for fat, lean, and bone mineral content.

This is a 2-equation, 3-unknown system. **It is mathematically underdetermined.**

### What Hologic does in practice

Hologic's APEX software applies a **soft-tissue constraint**: in pixels
classified as soft-tissue only (no bone), fat is assumed absent of BMC, reducing
the system to 2 unknowns. In bone pixels, an additional constraint is imposed
(typically: the pixel is modelled as a two-compartment mixture of soft tissue
and bone mineral, using pre-solved soft-tissue values from adjacent non-bone
pixels).

**Do not use a naive 2×3 pseudo-inverse** to attempt to invert the full
fat/lean/BMC system from two measurements alone — it will produce numerically
unstable, physically meaningless results. The correct approach requires the
region-prior constraints and calibrated coefficients that APEX uses internally.

---

## 5. Sensitivity to Repositioning

**Confirmed from DXA literature.**

DXA regional measurements — particularly for sub-regions of the body (e.g.
appendicular lean mass, arm fat fraction) — are sensitive to small differences
in:

- Patient position and limb rotation.
- Scan plane alignment.
- Inclusion or exclusion of border pixels at region boundaries.

Reported reproducibility coefficients of variation (CV) in the literature range
from approximately 1–3% for whole-body composition, but can exceed 5% for
specific sub-regions when repositioning protocols are not strictly followed.

**Consequence:** comparisons across sessions (e.g. longitudinal tracking with
BodyLoop) must account for repositioning variability. The pipeline does not
correct for this source of error.

---

## 6. 2D Projection: CoM and Inertia Are Planar, Not 3D

**Confirmed — inherent to DXA physics.**

DXA produces a 2D projection image: all tissue depth is collapsed along the
X-ray beam axis. Consequently:

- **Center of mass (CoM)** computed from a DXA map is a projected 2D CoM in
  the scan plane, not the true 3D CoM of the anatomical segment.
- **Moment of inertia** (planar moment) is the second moment of projected areal
  density about the beam-perpendicular axes. It does not equal the principal
  moment of inertia of the 3D body segment.

The true 3D inertia tensor has six independent components (three principal
moments, three products of inertia). DXA provides, at best, one planar moment
about each of the two scan-plane axes. The third axis (along the beam) cannot
be recovered from 2D projection data.

**Consequence:** mechanical models that require 3D inertia tensors (e.g. inverse
dynamics for BodyLoop gait analysis) must not use DXA-derived planar moments
directly without appropriate assumptions about the third-axis distribution.

---

## 7. Fan-Beam Geometry: Pixel Area Is Not Constant

**Confirmed — hardware design constraint.**

Hologic uses a fan-beam X-ray geometry. Unlike a pencil-beam scanner where all
pixels subtend equal solid angles, a fan-beam scanner produces pixels whose
effective area increases with distance from the beam isocenter. Peripheral pixels
cover a larger area of the patient than central pixels.

**Consequence:**

- Scalar pixel area approximations (computed from `PixelSpacing` × `PixelSpacing`)
  are valid only at the isocenter.
- Peripheral pixel masses are underestimated if the scalar approximation is used.
- The pipeline's `MapGeometry.pixel_area_cm2_scalar()` is provided for
  pencil-beam geometry only and is labelled accordingly.
- When calibrated fan-beam correction is available, use the per-pixel
  `pixel_area_cm2` array in `QuantitativeMap`.

---

## 8. Posture Differences: DXA (Supine) vs. BodyLoop (Standing)

**Confirmed — systematic experimental design issue.**

DXA scans are acquired with the patient lying supine on a flat table. BodyLoop
motion capture data are acquired with the participant standing. This difference
has several consequences:

- Soft-tissue gravitational pooling is different between supine and standing.
- Segment length in the scan axis changes slightly due to joint load bearing.
- The projected CoM in the coronal plane (DXA) does not equal the standing CoM
  in the same plane.

**There is no validated method to correct DXA-derived inertial parameters from
supine to standing geometry** without additional 3D imaging (CT or MRI).
Results from this pipeline should be interpreted as supine-geometry estimates.

---

## Summary Table

| Limitation | Status | Blocks pipeline? |
|---|---|---|
| SC prohibition | Confirmed | Yes (enforced) |
| Proprietary P/R format | Confirmed | Yes (cannot decode) |
| Missing calibration | Confirmed | Yes (quantitative reconstruction blocked) |
| Underdetermination (3 unknowns, 2 equations) | Confirmed | Yes (no naive pseudo-inverse) |
| Repositioning sensitivity | Confirmed | No (must be handled externally) |
| 2D projection (CoM/inertia) | Confirmed | No (results are 2D estimates) |
| Fan-beam pixel area variation | Confirmed | No (use per-pixel area when available) |
| Supine vs. standing posture | Confirmed | No (interpret as supine estimates) |
