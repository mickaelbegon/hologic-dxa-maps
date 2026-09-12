# Security Policy

## Scope

This project processes **medical imaging data** (DXA DICOM). Security issues include:

- Exposure of protected health information (PHI) through logs, outputs, or error messages
- Path traversal vulnerabilities when handling DICOM tag values as file names
- Malicious DICOM files triggering code execution
- Dependency vulnerabilities affecting data integrity

## Reporting a Vulnerability

Please **do not** open a public GitHub issue for security vulnerabilities.

Report to: mickael.begon@umontreal.ca  
Subject: `[SECURITY] hologic-dxa-maps — <brief description>`

Include:
- Description of the vulnerability
- Steps to reproduce (with synthetic/anonymized data only)
- Potential impact
- Suggested fix if any

Expected response time: within 5 business days.

## Data Handling Principles

- No patient data should ever be committed to this repository
- Logs must never contain patient identifiers
- File names derived from DICOM tags must be sanitized before use as filesystem paths
- The `.gitignore` is configured to block common medical data file extensions
