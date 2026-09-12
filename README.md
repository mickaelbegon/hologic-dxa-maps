# hologic-dxa-maps

[![CI](https://github.com/mickaelbegon/hologic-dxa-maps/actions/workflows/ci.yml/badge.svg)](https://github.com/mickaelbegon/hologic-dxa-maps/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)

Reproducible Python pipeline for extracting and reconstructing pixel-level quantitative maps from **Hologic DXA DICOM** files.

Developed for a project creating a bank of anthropometric and inertial segmental parameters in acrobats.

---

## The scientific problem

Dual-energy X-ray absorptiometry (DXA) measures the attenuation of X-rays at two energies to decompose body composition pixel by pixel into:

- **σ_fat(i,j)** — fat areal density [g/cm²]
- **σ_lean(i,j)** — lean tissue areal density [g/cm²]
- **σ_BMC(i,j)** — bone mineral content areal density [g/cm²]
- **σ_total(i,j) = σ_fat + σ_lean + σ_BMC**

From these maps one can compute:
- Regional and whole-body masses
- Projected centre of mass
- Planar moments of inertia (about the axis normal to the detector)

This pipeline enables the construction of **patient-specific segmental inertial parameters** for biomechanical modelling of acrobats.

---

## Why a Secondary Capture image is not enough

Hologic DXA scanners produce several DICOM object types. **Display-only images (Secondary Capture, SOPClassUID 1.2.840.10008.5.1.4.1.1.7) cannot be used to extract quantitative data.** Their pixel values are:

1. 8-bit or 16-bit display-scaled integers with no documented physical unit
2. Post-processed for visual appearance (contrast enhancement, log compression)
3. Not linked to any calibration reference

**Converting Secondary Capture pixel values to g/cm² is scientifically invalid.** This pipeline enforces this constraint and will raise an error if such a conversion is attempted.

---

## Hologic DICOM object types

| Object type | SOPClassUID (approx.) | Contains quantitative data? |
|---|---|---|
| Secondary Capture | 1.2.840.10008.5.1.4.1.1.7 | **No — display only** |
| Hologic Archive (DICOMDIR) | proprietary | Yes, in P/R files (proprietary format) |
| Structured Report (APEX results) | 1.2.840.10008.5.1.4.1.1.88.x | Yes — regional values with codes and units |
| Parametric Map | 1.2.840.10008.5.1.4.1.1.30 | Yes — pixel-level with RWVM |

### P and R files

Hologic archive DICOM objects embed two proprietary binary files (P and R) in private tags of group `(0023,xxxx)`. The structure of these files is **not publicly documented**. This pipeline can:

- Detect and extract the raw bytes (byte-for-byte, SHA-256 verified)
- Perform forensic inspection (size, entropy, magic bytes)
- **NOT decode** them without the Hologic SDK or official documentation

---

## Data capability levels

The software explicitly separates:

| Level | Description | Requires |
|---|---|---|
| **DICOM inspection** | Classify and inventory all objects | Any DICOM file |
| **P/R extraction** | Extract raw proprietary bytes | Archive DICOM |
| **SR extraction** | Parse regional numerical results | Structured Report DICOM |
| **Map import** | Load pre-calibrated arrays | Calibrated export (`.npy` + `metadata.json`) |
| **Experimental reconstruction** | Attempt bi-energy decomposition | P/R documentation + calibration phantoms |
| **BodyLoop fusion** | Register 3D surface with DXA projection | BodyLoop mesh + registration parameters |

---

## Installation

### Prerequisites

- Python 3.11 or newer
- pip ≥ 24

### From source (recommended for development)

```bash
git clone https://github.com/mickaelbegon/hologic-dxa-maps.git
cd hologic-dxa-maps
pip install -e ".[dev]"
pre-commit install
```

### With conda/miniforge

```bash
conda env create -f environment.yml
conda activate hologic-dxa
pip install -e ".[dev]"
pre-commit install
```

### Windows note

All dependencies install on Windows via pip without Fortran compilers. The optional `open3d` package (BodyLoop module) is not available on Windows via conda — use pip.

---

## Quick start

```bash
# Check pipeline status
hologic-dxa doctor

# Inspect a directory of DICOM files
hologic-dxa inspect /path/to/dicoms --output output/inventory

# Extract Hologic P/R files from archive DICOMs
hologic-dxa extract-pr /path/to/dicoms --output output/extracted

# Parse Structured Report results
hologic-dxa parse-sr /path/to/sr_dicoms --output output/sr

# Validate mass conservation (requires calibrated maps)
hologic-dxa validate-maps output/maps/bundle.h5 \
    --sr output/sr/sr_long_format.csv \
    --output output/validation
```

---

## Environment variables

| Variable | Description | Required |
|---|---|---|
| `DICOM_INPUT_DIR` | Root directory for DICOM search | No |
| `DICOM_SR_DIR` | Structured Report DICOM directory | No |
| `APEX_CSV_PATH` | APEX CSV export | No |
| `HOLOGIC_SDK_PATH` | Hologic SDK library path | No |
| `HOLOGIC_MAP_EXPORT_DIR` | Pre-calibrated map arrays | No |
| `BODYLOOP_POINT_CLOUD` | BodyLoop mesh or point cloud | No |
| `OUTPUT_DIR` | Output base directory | No (default: `output/`) |

---

## Privacy and data protection

- **No patient data is ever committed to this repository** (`.gitignore` blocks all DICOM extensions)
- Logs never contain patient names, IDs, birth dates, or accession numbers
- Research pseudonyms (`SUB_XXXXXXXX`) replace patient identifiers in all outputs
- DICOM source files are never modified
- De-identified copies are placed in a separate directory
- Private tag extraction must precede any de-identification that could remove P/R data

See [`docs/privacy.md`](docs/privacy.md) for full details.

---

## Synthetic reproducible example

```python
from tests.synthetic_data.generators import make_hologic_archive, make_parametric_map
import numpy as np

# Create a synthetic parametric map (no patient data)
ds = make_parametric_map(
    rows=50, cols=40,
    float_values=np.random.uniform(50, 200, (50, 40)).astype(np.float32),
    unit_code="mg/cm2",
    pixel_spacing_mm=(2.0, 2.0),
)
print(ds.SOPClassUID)  # 1.2.840.10008.5.1.4.1.1.30
```

See [`examples/synthetic_demo.ipynb`](examples/synthetic_demo.ipynb) for a full end-to-end example.

---

## Running tests

```bash
# Unit and property tests (no patient data required)
pytest tests/ -m "not real_data" -v

# With coverage
pytest tests/ -m "not real_data" --cov=hologic_dxa --cov-report=html

# Include slow tests
pytest tests/ -m "not real_data" --run-slow
```

---

## Limitations

See [`docs/limitations.md`](docs/limitations.md) for the complete list. Key constraints:

1. **P/R format is proprietary** — decoding requires Hologic SDK or documentation
2. **Two energies, three unknowns** — the system is underdetermined for bone pixels; no naive pseudo-inverse
3. **2D projection** — CoM and inertia are planar; 3D values require BodyLoop fusion
4. **Fan-beam geometry** — pixel area is not constant across the image
5. **Posture difference** — DXA is supine, BodyLoop is standing; registration is non-trivial

---

## Citation

If you use this software in your research, please cite it using [`CITATION.cff`](CITATION.cff).

---

## Contributing

Please read [`SECURITY.md`](.github/SECURITY.md) before reporting any vulnerability.

Bug reports: [open an issue](https://github.com/mickaelbegon/hologic-dxa-maps/issues/new/choose).
