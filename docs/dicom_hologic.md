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

## 7. P File Binary Format (Empirical — Hologic Horizon W, APEX 13.6)

**[CONFIRMED]** The following structure was determined by reverse engineering
a real Hologic Horizon W whole-body scan (APEX 13.6.0.7, serial 300829M).
Results may differ for other scanner models or APEX versions.

### 7.1 TLV Envelope

The P file uses a **TLV (Type-Length-Value) stream**:

```
  record := type(uint16 LE)  +  len(uint32 LE)  +  payload(len bytes)
```

The `len` field is the **raw payload size** (NOT including the 6-byte header).

### 7.2 Known P File Record Types

| Type | Hex | Size | Content |
|------|-----|------|---------|
| 1001 | 0x03E9 | ~26 B | Scanner serial number (UTF-16-LE, null-terminated) |
| 70   | 0x0046 | 4,259,907 B | Main scan data (sub-header + raw pixel counts) |

### 7.3 Record Type=70 Structure

```
  [6-byte TLV header]
  [67-byte sub-header]
    - bytes  0–33: Scan GUID (UTF-16-LE, 17 chars, e.g. "0396A55A49D74A69F")
    - bytes 34–66: zeros / padding
  [4,259,840 bytes of raw scan data]
    - shape: (16640, 128) uint16, row-major
    - lines 0, 2, 4, … → Low energy (100 kV)
    - lines 1, 3, 5, … → High energy (140 kV)
    → arr_L = data[0::2, :]  shape (8320, 128)
    → arr_H = data[1::2, :]  shape (8320, 128)
```

Size check: 67 + 16640 × 128 × 2 = 67 + 4,259,840 = **4,259,907** ✓

### 7.4 Blob Rest (after record type=70)

After the type=70 payload, the P file contains 7,467,413 more bytes:

```
  [3,477 bytes of zeros]  ← alignment / end-of-record padding
  [7,463,936 bytes of float64 data]
    - shape: (7289, 128), dtype=float64 (little-endian)
    - row 0: step-wedge attenuation curve (detectors 0-23 = 0.177→2.165,
             detectors 24-127 = 0.0)
    - subsequent rows: sparse attenuation map, predominantly 0 for air,
      increasing values for tissue and bone
```

**[SUSPECTED]** This float64 array represents the processed log-attenuation
after drum calibration, used internally by APEX. The 7289 rows may correspond
to the 8320 raw L-energy lines with ~1031 drum-band lines removed.

**[UNKNOWN]** The exact algorithm mapping (arr_L, arr_H, drum_ref) → float64
attenuation map. Requires Hologic SDK to implement correctly.

---

## 8. R File Binary Format (Empirical — Hologic Horizon W, APEX 13.6)

**[CONFIRMED]** The R file uses a **different TLV convention** from the P file:

```
  record := type(uint16 LE)  +  len_total(uint32 LE)  +  payload
  payload_size = len_total - 6   ← len_total INCLUDES the 6-byte header
```

### 8.1 Full R File Record Inventory (APEX 13.6, whole-body scan)

| Type | Hex | Payload | Content |
|------|-----|---------|---------|
| 90  | 0x005A | 6 B  | Version/flags: `01 05 02 01 01 01` |
| 223 | 0x00DF | 6 B  | Flags: `01 01 02 02 01 01` |
| 91  | 0x005B | 6 B  | Flags: `01 01 02 02 01 01` |
| 56  | 0x0038 | 11 B | ACF string: `"  1.024217"` (null-terminated) |
| 57  | 0x0039 | 11 B | Reference value: `" 13.030500"` |
| 58  | 0x003A | 2 B  | uint16 = 636 (2 × ROI Width) |
| 59  | 0x003B | 2 B  | uint16 = 150 (ROI Height — matches SR DICOM) |
| 60  | 0x003C | 2 B  | uint16 = 1 |
| 61  | 0x003D | 2 B  | uint16 = 4 |
| 62  | 0x003E | 2 B  | uint16 = 4 |
| 69  | 0x0045 | 2 B  | uint16 = 0 |
| 291 | 0x0123 | 7 B  | Protocol string: `"A237.6"` (null-terminated) |
| 221 | 0x00DD | 144 B| Structured header (indices, offsets) |
| 202 | 0x00CA | 190,800 B | **95,400 uint16** — calibration data (150×636) |
| 224 | 0x00E0 | 190,800 B | **95,400 uint16** — calibration data (150×636) |
| 225 | 0x00E1 | 11 B | Calibration value: `"  1.921908"` |
| 251 | 0x00FB | 2 B  | uint16 = 225 |
| 250 | 0x00FA | 2 B  | uint16 = 6 |
| 126 | 0x007E | 60 B | Structured descriptor (offsets, sizes) |
| 253 | 0x00FD | 8 B  | float64 = 0.0902 |
| 254 | 0x00FE | 8 B  | float64 = 0.8687 |
| 252 | 0x00FC | 144 B| Structured header (indices, offsets) |
| 222 | 0x00DE | 2,457,000 B | **Large calibration / segmentation block** |
| 0   | 0x0000 | 0 B  | End-of-file sentinel |

### 8.2 Notes on Key Records

**Records 202 and 224** (190,800 bytes each): Contain 95,400 uint16 values
each, with virtually identical distributions (ratio type202/type224 ≈ 1.000).
**[SUSPECTED]** These are L and H drum-calibration reference counts sampled at
each scan position, stored as (150, 636) arrays — **not** the bi-energy body
scan images. Values range 221–2368, median ≈ 600.

**Record 222** (2,457,000 bytes): Large block whose first bytes contain uint16
values in the 220–225 range (matching R file type codes), followed by larger
values. **[UNKNOWN]** Purpose — possibly an encoded segmentation map, a
compressed pixel classification array, or a nested record structure.

**Records 253 and 254** (float64): Physical calibration constants
(0.0902 and 0.8687). **[SUSPECTED]** These are the BCF (bone calibration
factor) intermediate values used by APEX for areal density computation.

**Calibration string records (56, 57, 225)**:
- type=56: `ACF = 1.024217` (Air Calibration Factor)
- type=57: `13.030500` (reference thickness or slope — matches APEX XML)
- type=225: `1.921908` (reference constant)

---

## 9. Scan Identification (GUID)

**[CONFIRMED]** Each Hologic scan is assigned a unique GUID that appears in
multiple locations:

| Location | Format | Example |
|----------|--------|---------|
| P file sub-header bytes 0–33 | UTF-16-LE, 17 chars | `"0396A55A49D74A69F"` |
| SR DICOM PatientName | `"FCA{GUID}^{hash}"` | `"FCA0396A55A49D74A69F^281790962FD450592773"` |

The `FCA` prefix in the SR PatientName identifies Hologic's anonymization
scheme. The 17-character hex GUID uniquely links the SR DICOM to the P file
embedded in the archive DICOM.

---

## 10. APEX XML Tag (0019,1000)

**[CONFIRMED]** Tag `(0019,1000)` contains a UTF-8 XML document with the full
APEX body composition table. Key fields observed:

```xml
<ScanID>A10152008</ScanID>
<k>1.181147</k>
<d0>51.406250</d0>
<ACF>1.023</ACF>
<BCF>0.987</BCF>
```

Plus per-region results: Area, BMC, BMD, Fat, Lean, %Fat for all body regions.
These values match the SR DICOM measurements exactly.

---

## 11. Open Questions Requiring Hologic Documentation

The following questions cannot be answered from empirical observation alone and
require either the Hologic SDK or official Hologic technical documentation to
resolve:

1. **P file blob_rest → density map**: The exact algorithm that maps
   (arr_L, arr_H, drum_ref) to the float64 attenuation array in blob_rest.
   Specifically: how are the ~1031 drum-band lines identified and removed,
   and what I₀ reference is used for `log(I₀/I)`.

2. **R file record 222**: The 2.46 MB block — is it a segmentation label
   image, a nested sub-record stream, or compressed body composition data?

3. **R file records 202/224**: Confirmed as nearly-identical uint16 arrays.
   Are they L and H drum reference counts, or some other calibration quantity?

4. **Encoding scheme versions**: What does each `ENCODING_SCHEME` string value
   imply about the binary format of P/R files?

5. **Multi-creator blocks**: Can group `(0023)` contain more than one private
   creator, and if so, how do the data tag offsets change?

6. **Scan-type tagging**: Is there a private tag that identifies the scan type
   (e.g. whole-body, spine, hip), beyond what can be inferred from
   `BodyPartExamined` or the series description?

7. **APEX 5.x compatibility**: Are the private tag definitions for APEX 5.x
   backward-compatible with the schema documented here?
