# Hologic DXA DICOM Specifics

This document describes what is known, suspected, and unknown about the DICOM
files produced by Hologic APEX software. All claims are based on empirical
observation of real DICOM files and publicly available documentation fragments.
Nothing in this document should be interpreted as official Hologic technical
documentation.

**Label convention throughout this document:**
- **[CONFIRMED]** — Verified against real DICOM files or official sources.
- **[SUSPECTED]** — Consistent with observed data but not independently verified.
- **[UNKNOWN]** — Cannot be determined from available information alone.

---

## 1. Private Tag Group (0023)

### 1.1 Private Creator Tag

| Tag | Value |
|-----|-------|
| `(0023,0010)` | `"HOLOGIC, Inc."` |

**[CONFIRMED]** The private creator string `"HOLOGIC, Inc."` is registered at
element `(0023,0010)`. All subsequent tags in this group are within the private
block reserved by this creator. The full tag addresses listed below assume
the creator is at offset `0x10` (the first available block), yielding data
element offsets in the range `(0023,1000)` through `(0023,10FF)`.

### 1.2 Known Tags Within Block (0023,10xx)

| Tag | Offset | VR (observed) | Name | Status |
|-----|--------|---------------|------|--------|
| `(0023,1000)` | `0x00` | `LO` | Encoding Scheme | **[SUSPECTED]** |
| `(0023,1001)` | `0x01` | `LO` | P File Name | **[SUSPECTED]** |
| `(0023,1002)` | `0x02` | `OB` | P File Data | **[CONFIRMED]** |
| `(0023,1003)` | `0x03` | `UL` | P File Length | **[CONFIRMED]** |
| `(0023,1004)` | `0x04` | `OB` | R File Data | **[CONFIRMED]** |
| `(0023,1005)` | `0x05` | `UL` | R File Length | **[CONFIRMED]** |

The presence of these tags depends on the APEX version and the type of scan
(whole-body, spine, hip, forearm). Not all tags appear in all archive objects.

### 1.3 Encoding Scheme Tag (0023,1000)

**[SUSPECTED]** This tag appears to contain a short string identifying the
encoding format used for the embedded P and R binary files. Observed values
include strings of the form `"HOLOGIC_ENCODING_Vx"`. The exact mapping from
scheme string to binary format is **[UNKNOWN]** — Hologic has not published
documentation for the internal encoding.

**Do not attempt to parse P/R file content based on this tag alone** without
reference to official Hologic SDK documentation.

---

## 2. Archive Objects vs. Display Objects

### 2.1 Identification Criteria

Hologic DXA studies typically contain two categories of DICOM objects:

**Archive (data) objects:**
- Contain private tags in group `(0023)` with non-empty P/R blobs.
- `PRDiagnosis` → `PR_PRESENT`.
- `DicomClass` → `HOLOGIC_ARCHIVE`.
- May have a Secondary Capture SOP Class UID while still being classified as an
  archive (the PR data takes precedence in classification).

**Display objects:**
- Use the Secondary Capture SOP Class (`1.2.840.10008.5.1.4.1.1.7`).
- Contain no private group `(0023)` data, or group `(0023)` is absent.
- Pixel values are windowed/scaled for display and are **not** physical units.
- `DicomClass` → `SECONDARY_CAPTURE`.

**Programmatic check:**

```python
from hologic_dxa.dicom.classify import classify, DicomClass

result = classify(ds)
if result.dicom_class == DicomClass.HOLOGIC_ARCHIVE:
    # Has P/R data; proceed to extraction (requires SDK or calibration)
    ...
elif result.dicom_class == DicomClass.SECONDARY_CAPTURE:
    # Display only — must not be used as a quantitative map
    ...
```

### 2.2 Distinction Is Not Always Clear

**[SUSPECTED]** Some APEX versions produce objects that have a Secondary Capture
SOP Class UID but also include private group `(0023)` tags. In these cases, the
pipeline uses the PR diagnosis to override the SOP-class-based classification
and marks the object as `HOLOGIC_ARCHIVE`.

**[UNKNOWN]** Whether this dual-encoding occurs only in specific APEX versions
or in specific scan modes (e.g. follow-up scans stored alongside the original).

---

## 3. APEX Version Variations

**[SUSPECTED]** The layout and completeness of private tags varies across APEX
software versions. Observations:

| APEX generation | P/R Tags | Notes |
|-----------------|----------|-------|
| APEX 3.x | Present in archive objects | Typical layout as documented above |
| APEX 4.x | Present; encoding scheme tag more frequently populated | — |
| APEX 5.x | **[UNKNOWN]** — insufficient samples | — |

The `SoftwareVersions` DICOM attribute (tag `0018,1020`) may identify the APEX
version string. This field is not always populated and its exact format varies.

**[UNKNOWN]** Whether private tag offsets shift between APEX versions if
Hologic registers multiple private creators in group `(0023)`. If a second
creator is registered at offset `0x11`, the data tags would move to
`(0023,1100)` through `(0023,11FF)`. The pipeline assumes a single creator at
offset `0x10`.

---

## 4. Modality and Manufacturer Fields

**[CONFIRMED]** Hologic DXA systems produce files with:

- `Manufacturer` = `"HOLOGIC, Inc."` (exact string, including the period)
- `Modality` = `"BMD"` for bone density scans or `"OT"` for Other

Some older APEX versions may use `"OT"` for whole-body composition scans that
are not primarily bone density acquisitions. The pipeline accepts `BMD`, `OT`,
and `DX` as valid DXA modalities.

---

## 5. Parametric Map Objects

Some APEX versions or processing pipelines produce DICOM Parametric Map objects
(SOP Class `1.2.840.10008.5.1.4.1.1.30`) with:

- `FloatPixelData` or `DoubleFloatPixelData` — float-valued pixel arrays.
- `RealWorldValueMappingSequence` — slope/intercept mapping to physical units.

**[CONFIRMED]** These objects are classified as `PARAMETRIC_MAP` and are
`ELIGIBLE` for quantitative use when both RWVM and float pixel data are present.

**[SUSPECTED]** The physical unit stored in `MeasurementUnitsCodeSequence`
within the RWVM sequence is one of `"mg/cm2"`, `"g/cm2"`, or `"HU"` (for
reconstructed CT-like values), depending on what the post-processing step
outputs.

**[UNKNOWN]** Whether parametric maps produced by third-party post-processing
tools (rather than APEX itself) carry the same private tag layout.

---

## 6. Structured Report Objects

Hologic APEX can export regional composition results as DICOM Structured Reports
(SR), typically using the Comprehensive SR SOP Class
(`1.2.840.10008.5.1.4.1.1.88.33`). These SRs contain:

- `ContentSequence` with `NUM` value type items.
- Each item carries a concept code (e.g. `"Bone Mineral Density"`), a numeric
  value, and a unit code.

**[CONFIRMED]** SR objects are classified as `STRUCTURED_REPORT` with
`QuantitativeEligibility.ELIGIBLE` because they may contain calibrated regional
measurements already computed by APEX.

**[SUSPECTED]** The coding scheme used is a mix of DICOM-defined codes (`DCM`)
and Hologic-specific local codes whose meanings are not publicly documented.

---

## 7. Open Questions Requiring Hologic Documentation

The following questions cannot be answered from empirical observation alone and
require either the Hologic SDK or official Hologic technical documentation to
resolve:

1. **P file internal format**: What is the byte-level structure of the P file?
   Is it a proprietary binary, a compressed format, or a structured record?

2. **R file purpose**: Is the R file a reference scan, a calibration record, or
   something else? Does it correspond to the bone reference phantom subtraction?

3. **Encoding scheme versions**: What does each `ENCODING_SCHEME` string value
   imply about the binary format of P/R files?

4. **Calibration coupling**: Are the calibration factors needed to convert P
   data to areal density stored in the R file, in a separate DICOM object, or
   in the scanner's internal database?

5. **Multi-creator blocks**: Can group `(0023)` contain more than one private
   creator, and if so, how do the data tag offsets change?

6. **Scan-type tagging**: Is there a private tag that identifies the scan type
   (e.g. whole-body, spine, hip), beyond what can be inferred from
   `BodyPartExamined` or the series description?

7. **APEX 5.x compatibility**: Are the private tag definitions for APEX 5.x
   backward-compatible with the schema documented here?
