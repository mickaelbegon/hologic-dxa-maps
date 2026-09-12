# Implementation Plan

## Phase dependencies

```
Phase 0: Infrastructure (repo, CI, config)
    │
    ├── Phase 1: DICOM inventory & classification    ← no dependencies
    ├── Phase 3: SR parsing                          ← no dependencies
    └── Phase 4: Data models                         ← no dependencies
         │
         ├── Phase 2: P/R extraction (needs Phase 1 + Phase 4)
         │
         └── Phase 5: Provider interface (needs Phase 4)
              │
              ├── Phase 6: Bi-energy architecture (needs Phase 5, BLOCKED on Hologic docs)
              │
              └── Phase 7: Quantitative calculations (needs Phase 5)
                   │
                   ├── Phase 8: Validation (needs Phase 7)
                   └── Phase 9: Export + CLI (needs Phase 5, 7, 8)
                        │
                        └── Phase 10: BodyLoop fusion (needs Phase 9, experimental)
```

## Current status

| Phase | Status | Blocking factors |
|---|---|---|
| 0 Infrastructure | ✅ Complete | — |
| 1 DICOM inventory | ✅ Core implemented | None |
| 2 P/R extraction | ✅ Implemented | Decoding still blocked (no SDK) |
| 3 SR parsing | 🔄 In progress | None |
| 4 Data models | ✅ Implemented | None |
| 5 Provider interface | 🔄 In progress | None |
| 6 Bi-energy reconstruction | 🔲 Architecture only | Hologic SDK / documentation |
| 7 Quantitative calculations | 🔄 In progress | None |
| 8 Validation | 🔄 In progress | Real calibrated maps |
| 9 Export + CLI | 🔄 In progress | None |
| 10 BodyLoop | 🔲 Interface stub | BodyLoop data + registration research |

## Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Hologic never provides SDK/docs for P/R | Medium | SR + APEX CSV provide regional values; exported maps from APEX can be imported |
| APEX SR format varies by version | High | Concept mapping YAML is user-editable; raw codes always preserved |
| Fan-beam geometry parameters not in DICOM | Medium | Document limitation; accept scalar pixel area with warning |
| BodyLoop ↔ DXA registration quality insufficient | Medium | Articulated + deformable registration; uncertainty quantification |
| Pixel area not constant (fan-beam) | Certain | `pixel_area_cm2` field accepts `np.ndarray`; document clearly |

## Information still needed from Hologic

To unblock Phase 6 (experimental bi-energy reconstruction):

1. P file binary format specification:
   - Native array dimensions (rows × cols × energy levels)
   - Numeric dtype (int16, uint16, float32, …)
   - Byte order (little/big endian)
   - Layout (interleaved energies, or separate planes)

2. R file binary format specification (same questions)

3. Dark current reference values or acquisition protocol

4. Air reference image acquisition protocol

5. Detector gain correction map format

6. Calibration function F_θ(σ_fat, σ_lean) → (a_L, a_H) coefficients or lookup table

7. For bone pixels: additional constraint or algorithm for the underdetermined system

8. Final unit conversion factors (ADU → g/cm²)

Any of the following would also unblock the phase:
- Access to the Hologic SDK with Python bindings
- Official export from Hologic software of calibrated pixel maps in an open format
- Calibration phantom with known composition + matching archive DICOM

## Testing strategy

- **No real patient data in the repository** — ever
- Synthetic DICOM generators cover all object types
- Property-based tests (hypothesis) cover mathematical invariants
- Integration tests with real data use `@pytest.mark.real_data` and are excluded from CI

## Milestone: minimal viable pipeline

The pipeline is considered functional when:
1. ✅ `hologic-dxa doctor` reports correctly
2. ✅ `hologic-dxa inspect` produces PHI-free inventory
3. ✅ Archive DICOMs are correctly classified and extracted
4. ✅ SR results are parsed with codes and units
5. ✅ Pre-calibrated map arrays can be imported and validated
6. ✅ Mass, CoM, and inertia are computed from a synthetic bundle
7. ✅ All tests pass with `pytest -m "not real_data"`
8. ✅ `ruff check` and `mypy` report no blocking errors
