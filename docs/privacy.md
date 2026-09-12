# Privacy and Data Protection

## Principles

This pipeline processes medical imaging data (DXA DICOM). The following principles are non-negotiable:

1. **No PHI in logs** — Patient names, IDs, birth dates, sex, weight, size, accession numbers, institution names, and physician names are never written to log files or output reports.
2. **No PHI in outputs** — Inventory CSV/JSON, validation reports, and HDF5 exports contain only pseudonymised research identifiers.
3. **No modification of source DICOMs** — Source files are always opened read-only. Extracted bytes are never modified.
4. **Pseudonymisation precedes de-identification** — Private tag extraction (P/R files) must occur before any de-identification step that might remove those tags.
5. **Separate de-identified copies** — De-identified DICOMs are written to a distinct output directory, never alongside originals.

## PseudonymRegistry

The `PseudonymRegistry` class (`hologic_dxa.provenance`) maintains a persistent, bijective mapping:

```
PatientID (PHI) → research_id (e.g. SUB_A3F7B291)
```

- Stored in a JSON file at a configurable path (default: `output/pseudonym_registry.json`)
- Never logged or printed
- Research IDs are randomly generated UUIDs (8 hex chars, uppercase prefix `SUB_`)
- The mapping is bijective: a collision raises `ValueError`
- The registry file must be stored outside the repository and protected with appropriate access controls

## DICOM de-identification standard

When de-identification is needed, follow **DICOM PS 3.15 Annex E** (Basic Application Level Confidentiality Profile). Key points:

- Remove or empty all Type 1 and Type 2 patient attributes
- Replace UIDs consistently (same original UID → same replacement UID within a study)
- Preserve private tags containing P/R data **before** applying de-identification

## What never appears in any output file

- PatientName
- PatientID
- PatientBirthDate
- PatientSex, PatientAge, PatientWeight, PatientSize
- OtherPatientIDs, OtherPatientNames
- InstitutionName, InstitutionAddress
- ReferringPhysicianName, PerformingPhysicianName
- RequestingPhysician, OperatorsName
- AccessionNumber

## Responsibility

Researchers using this pipeline are responsible for:
- Securing the `pseudonym_registry.json` file
- Ensuring output directories are access-controlled
- Complying with their institutional ethics approval and applicable data protection regulations (PIPEDA, GDPR, etc.)
