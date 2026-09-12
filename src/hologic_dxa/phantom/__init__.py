"""Phantom calibration support.

This module bridges the gap between 'experimental reconstruction' (Phase 6)
and 'validated calibration'. Phantom acquisitions with known composition
provide the ground truth needed to fit calibration functions F_θ.

Currently implemented:
- Phantom data loading interface
- Expected reference values for common Hologic phantoms

Not yet implemented (requires phantom data):
- Calibration function fitting
- Uncertainty propagation from calibration to maps
"""
