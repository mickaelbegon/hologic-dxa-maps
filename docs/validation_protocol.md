# Validation Protocol

## Purpose

Validate that pixel-integrated masses from quantitative maps agree with APEX regional values from Structured Reports or CSV exports.

## Mass conservation check

For each tissue component and each anatomic region:

```
M_integrated = Σ_{i,j ∈ region} σ(i,j) × A_{ij}
```

Compared to the APEX reference `M_APEX` from the Structured Report.

Metrics reported:
- Absolute error: `|M_integrated - M_APEX|` [g]
- Relative error: `(M_integrated - M_APEX) / M_APEX` [dimensionless]
- Bland-Altman mean: `(M_integrated + M_APEX) / 2`

### Acceptance thresholds

**The thresholds below are PROVISIONAL.** They are placeholders and must be replaced by empirically derived values from:
- Repeated acquisitions on calibration phantoms (spine phantom, body composition phantom)
- Test-retest data from the same subjects
- Comparison across APEX software versions

| Metric | Provisional threshold |
|---|---|
| Relative error (whole body total mass) | ≤ 5% |
| Absolute error (whole body total mass) | ≤ 50 g |
| Relative error (regional) | ≤ 10% |

These thresholds are configurable in `config.py` and in the CLI `--mass-tolerance-relative` flag.

## Regional comparison

Regions expected from Hologic APEX Structured Reports:
- Whole body
- Head
- Trunk
- Left arm / Right arm
- Left leg / Right leg
- Sub-regions (pelvis, ribs, spine, etc.) when available

For each region, compare:
- Fat mass [g]
- Lean mass [g]
- BMC [g]
- Total mass [g]

## Phantom calibration

To establish thresholds and validate the pipeline, the following phantoms are recommended:
- **Hologic spine phantom** (included with every clinical Hologic system) — for BMD calibration
- **Hologic body composition phantom** — for soft tissue calibration
- **Anthropomorphic phantom** (e.g. CIRS Model 077) — for whole-body composition

Phantom acquisition protocol:
1. Acquire phantom on the same device used for subject scans
2. Run the pipeline on phantom DICOMs
3. Compare integrated masses to certified reference values
4. Repeat across scan positions (±2 cm repositioning)

## Repeatability assessment

For test-retest reliability:
1. Scan the same phantom (or volunteer) twice with repositioning
2. Compute integrated masses for both scans
3. Report coefficient of variation (CV%) and repeatability coefficient (RC)

## Output files

```
output/validation/
├── validation_results.csv   # per-region, per-component comparisons
├── validation_summary.json  # pass/fail, thresholds used, provenance
└── validation_report.html   # maps, residuals, Bland-Altman plots
```

The HTML report includes:
- Fat, lean, BMC maps with colour scale
- Residual maps (integrated - reference)
- Horizontal and vertical pixel profiles
- Regional comparison table
- All warnings and provenance metadata

## What is NOT validated here

- 3D inertia tensor (not directly measurable from DXA)
- Centre of mass in depth (Z-axis)
- Absolute accuracy of areal density values (requires phantom with certified composition)
- Agreement between APEX versions (requires multi-version dataset)
