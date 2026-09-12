# Architecture

## Overview

The pipeline is structured as a set of independent, composable layers. Each layer has explicit input/output contracts and fails explicitly when required data is absent.

```
┌─────────────────────────────────────────────────────────────────┐
│  CLI (hologic_dxa.cli)  — typer + rich                          │
└───────────────────────────────┬─────────────────────────────────┘
                                │
        ┌───────────────────────┼────────────────────────┐
        │                       │                        │
┌───────▼──────────┐  ┌─────────▼────────┐  ┌──────────▼────────┐
│ DICOM layer      │  │ Provider layer   │  │ Validation layer  │
│ dicom/inventory  │  │ providers/base   │  │ validation/       │
│ dicom/classify   │  │ providers/...    │  │ mass_conservation │
│ dicom/private_   │  │ (Protocol)       │  │ regional_compare  │
│ tags             │  │                  │  │ report            │
│ dicom/extract_pr │  └─────────┬────────┘  └──────────┬────────┘
│ dicom/sr         │            │                       │
└──────────────────┘            │                       │
                                ▼                       │
                    ┌───────────────────────┐           │
                    │ Maps layer            │           │
                    │ maps/bundle           │◄──────────┘
                    │ maps/geometry         │
                    │ maps/integration      │
                    │ maps/units            │
                    │ maps/quality_control  │
                    └───────────┬───────────┘
                                │
                    ┌───────────▼───────────┐
                    │ Export layer          │
                    │ export/hdf5           │
                    │ export/zarr           │
                    │ export/csv            │
                    │ export/dicom_pm       │
                    └───────────────────────┘
```

## Provider pattern

All sources of quantitative maps implement the `QuantitativeMapProvider` Protocol:

```python
class QuantitativeMapProvider(Protocol):
    name: str
    experimental: bool
    def can_read(self, source: Path) -> bool: ...
    def load(self, source: Path) -> QuantitativeMapBundle: ...
    def describe_capabilities(self) -> dict: ...
```

Adding a new data source (e.g. a new APEX version exporter) requires only implementing this protocol and registering the provider — no changes to the rest of the pipeline.

## Failure modes

The pipeline uses explicit failure rather than silent degradation:

| Missing data | Result |
|---|---|
| Pixel area / pixel spacing | `ValueError` — mass calculation refuses to run |
| Calibration factors | `QuantitativeDataUnavailableError` |
| P/R decoding parameters | `NotImplementedError` with actionable message |
| Hologic SDK | `NotImplementedError` with actionable message |
| SR with no numeric items | Empty `list[SRMeasurement]` + warning log |
| Mismatched map shapes | `ValueError` from `QuantitativeMapBundle.__post_init__` |

## Provenance chain

Every output carries:
1. SHA-256 of all input files
2. SOP Instance UID, Study UID, Series UID
3. Device manufacturer, model, APEX version
4. Pipeline version and timestamp
5. Ordered list of processing operations

This allows tracing any result back to its exact source files.

## Data flow for quantitative maps

```
DICOM files
    │
    ▼
[inspect] → inventory.json / inventory.csv
    │
    ├── Archive DICOMs
    │       │
    │       ▼
    │   [extract-pr] → {research_id}_{uid}_P.bin + manifest.json
    │
    ├── Structured Report DICOMs
    │       │
    │       ▼
    │   [parse-sr] → sr_long_format.csv + sr_tree.json
    │
    └── Parametric Map DICOMs  ──┐
                                 │
Pre-calibrated arrays ───────────┤
(exported_arrays provider)       │
                                 ▼
                    QuantitativeMapBundle (validated)
                                 │
                    ┌────────────┼────────────┐
                    ▼            ▼            ▼
               [validate]   [compute]    [export-hdf5]
               mass cons.   mass props   output.h5
```
