# Quantitative Map Framework

This document describes the physical quantities, mathematical formulas, and
geometric assumptions used in the `hologic-dxa-maps` pipeline when computing
derived quantities from DXA areal density maps.

---

## 1. Physical Quantities: σ_fat, σ_lean, σ_BMC

The pipeline works with three **areal density** maps, each denoted $\sigma$
(lower-case sigma) with a tissue subscript:

| Symbol | Quantity | Unit | Physical meaning |
|--------|----------|------|-----------------|
| $\sigma_{\text{fat}}$ | Fat areal density | g/cm² | Projected mass of fat tissue per unit scan-plane area |
| $\sigma_{\text{lean}}$ | Lean areal density | g/cm² | Projected mass of lean (non-fat, non-bone) tissue per unit scan-plane area |
| $\sigma_{\text{BMC}}$ | Bone mineral content areal density | g/cm² | Projected mass of bone mineral per unit scan-plane area |

All three quantities are **2D projections**: the scanner X-ray beam passes
through the entire depth of tissue, and what is measured is the depth-integrated
mass per unit area in the scan plane (not volumetric density in g/cm³).

### Relationship to total areal density

$$\sigma_{\text{total}}(x, y) = \sigma_{\text{fat}}(x, y) + \sigma_{\text{lean}}(x, y) + \sigma_{\text{BMC}}(x, y)$$

This additive decomposition holds when the three components are assumed to be
spatially co-located but non-overlapping in depth — a simplifying model used
by Hologic APEX. At bone pixels the model is more nuanced because bone and
soft tissue genuinely coexist along the beam path (see `docs/limitations.md`,
section 4).

---

## 2. Why pixel_area_cm² Is Required for Mass Integration

A pixel's **mass** is not its areal density value alone. It is:

$$m_i = \sigma_i \cdot A_i$$

where $\sigma_i$ is the areal density at pixel $i$ (g/cm²) and $A_i$ is the
effective area that pixel $i$ covers in the scan plane (cm²). Without $A_i$,
the sum $\sum_i \sigma_i$ has units of g/cm² (a mean, not a total mass) and is
physically meaningless as an absolute mass estimate.

The `QuantitativeMap` class therefore requires `pixel_area_cm2` to compute
`integrated_mass_g()` and raises `ValueError` if it is `None`.

---

## 3. Fan-Beam Geometry and Non-Uniform Pixel Area

Hologic DXA systems use **fan-beam** X-ray geometry: the X-ray source emits a
fan-shaped beam that subtends a fixed angular range. The detector is a line
array perpendicular to the beam.

Because the beam diverges, a pixel at the edge of the field of view corresponds
to a larger solid angle than a central pixel. The effective area $A_i$ is
**not constant** across the image:

$$A_i \approx \frac{\Delta\theta_i \cdot \Delta y}{\cos^2\theta_i} \cdot d^2$$

where $\Delta\theta_i$ is the angular width of pixel $i$, $d$ is the source-to-
isocenter distance, and $\theta_i$ is the angle of pixel $i$ from the central
beam. For off-center pixels at typical DXA field widths (~60 cm), the area
deviation from the central-pixel value can reach 5–15%.

**Consequence for the pipeline:**

- `MapGeometry.pixel_area_cm2_scalar()` returns a single scalar area computed
  from `PixelSpacing`. This is valid **only at the isocenter** and is labelled
  as a pencil-beam approximation.
- When fan-beam-corrected per-pixel areas are available (from the calibration
  step), they are stored as a 2D array in `QuantitativeMap.pixel_area_cm2`.
  Always prefer the array form over the scalar approximation.

---

## 4. Mathematical Formulas

### 4.1 Total Regional Mass

$$M = \sum_{i \in \mathcal{R}} \sigma_i \cdot A_i$$

where $\mathcal{R}$ is the set of pixels in the anatomical region of interest
(defined by `valid_mask` or an explicit mask passed to `regional_mass_g`).

Units: g (grams), since $[\text{g/cm}^2] \times [\text{cm}^2] = [\text{g}]$.

### 4.2 Center of Mass (Projected, 2D)

The 2D projected center of mass in the scan plane is:

$$\bar{x} = \frac{\sum_{i \in \mathcal{R}} \sigma_i \cdot A_i \cdot x_i}{\sum_{i \in \mathcal{R}} \sigma_i \cdot A_i}$$

$$\bar{y} = \frac{\sum_{i \in \mathcal{R}} \sigma_i \cdot A_i \cdot y_i}{\sum_{i \in \mathcal{R}} \sigma_i \cdot A_i}$$

where $(x_i, y_i)$ are the scan-plane coordinates of pixel $i$ centre in mm,
computed from the pixel index and `PixelSpacing`:

$$x_i = \text{col}_i \times \Delta x_{\text{mm}}, \quad y_i = \text{row}_i \times \Delta y_{\text{mm}}$$

Units: mm. This is the **projected** CoM — see section 5 for the distinction
from the 3D CoM.

In `center_of_mass_mm`, pixel spacing is provided as `pixel_spacing_mm = (row_spacing, col_spacing)`
and the origin is taken at pixel `(0, 0)`.

### 4.3 Planar Moment of Inertia

The **planar moment of inertia** (second moment of projected areal density about
the CoM) is:

$$I_{\perp} = \sum_{i \in \mathcal{R}} \sigma_i \cdot A_i \cdot r_i^2$$

where $r_i$ is the distance from pixel $i$ to the projected CoM:

$$r_i^2 = (x_i - \bar{x})^2 + (y_i - \bar{y})^2$$

Units: g·mm². This is analogous to the mass moment of inertia about an axis
perpendicular to the scan plane and passing through the projected CoM.

> **Note:** $I_{\perp}$ is NOT the same as the principal moment of inertia
> $I_{zz}$ from a 3D inertia tensor. See section 5.

---

## 5. 2D Projected Measurements vs. 3D Inertia Tensor

The 3D inertia tensor of a rigid body has six independent components:

$$\mathbf{I} = \begin{pmatrix} I_{xx} & I_{xy} & I_{xz} \\ I_{xy} & I_{yy} & I_{yz} \\ I_{xz} & I_{yz} & I_{zz} \end{pmatrix}$$

where $I_{xx} = \sum m_i (y_i^2 + z_i^2)$, etc.

DXA integrates along the beam axis (say, $z$). The projected areal density
$\sigma(x,y) = \int \rho(x,y,z)\,dz$ collapses the $z$-dimension. From DXA
alone:

- $I_{zz}^{\text{DXA}} \approx \sum \sigma_i A_i (x_i^2 + y_i^2)$ is
  computable (the planar moment $I_\perp$ above, shifted by the parallel-axis
  theorem from the segment endpoint rather than the CoM).
- $I_{xx}$, $I_{yy}$, $I_{xz}$, $I_{yz}$, $I_{xy}$ require knowledge of the
  3D depth distribution $\rho(x,y,z)$, which is lost in the projection.

**Consequence:** multi-body dynamics models requiring the full inertia tensor
must supplement DXA-derived quantities with assumptions about the 3D depth
distribution (e.g. cylindrical or ellipsoidal segment cross-sections). The
pipeline does not automate this step.

---

## 6. Calibrated vs. Display-Scaled Maps

### Calibrated maps

A **calibrated** map has pixel values in physical units (g/cm²) traceable to a
known bone-equivalent phantom. These are produced by:

1. Applying beam-hardening correction to raw detector signal.
2. Solving the two-energy system for each pixel (with appropriate tissue prior).
3. Applying phantom calibration factors to convert from equivalent water/bone
   thickness to areal density.

Calibrated maps are identified in DICOM by:
- A `RealWorldValueMappingSequence` with `MeasurementUnitsCodeSequence` citing
  a physical unit (`"g/cm2"`, `"mg/cm2"`).
- Or by being produced by a verified calibration step in this pipeline.

### Display-scaled maps

A **display-scaled** map has pixel values transformed for visual presentation:
- Windowed to the display bit depth.
- Level-shifted (rescale slope/intercept applied for display, not for physics).
- Possibly colour-mapped.

Display-scaled maps are **not** in physical units. Applying mass integration
formulas to display pixel values produces numbers with no physical meaning.

**Secondary Capture objects are always display-scaled** and are rejected by the
pipeline at the classification step (`QuantitativeEligibility.NOT_ELIGIBLE`).

---

## 7. Summary of Pipeline Assumptions

| Assumption | Consequence if violated |
|------------|------------------------|
| Pixel values are calibrated areal density (g/cm²) | Mass, CoM, and inertia will be dimensionally incorrect |
| `pixel_area_cm2` is supplied | `ValueError` from `integrated_mass_g()` and `regional_mass_g()` |
| Uniform pixel area (scalar) is valid | Peripheral pixel masses will be underestimated (fan-beam error) |
| Maps are co-registered (same grid) | Bundle validation fails with shape mismatch |
| Patient supine, anatomical position | CoM and inertia are supine-geometry estimates only |
