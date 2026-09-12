"""Shared pytest fixtures for the hologic-dxa test suite.

All fixtures produce synthetic datasets with no PHI.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.synthetic_data.generators import (
    make_empty_pr_archive,
    make_hologic_archive,
    make_parametric_map,
    make_secondary_capture,
    make_structured_report,
    make_truncated_pr_archive,
    save_to_tmpdir,
)


# ---------------------------------------------------------------------------
# In-memory DICOM dataset fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def secondary_capture_ds():
    """Minimal Secondary Capture dataset with synthetic pixel data."""
    return make_secondary_capture()


@pytest.fixture
def hologic_archive_ds():
    """Archive dataset containing both P and R file blobs."""
    return make_hologic_archive()


@pytest.fixture
def hologic_archive_padded_ds():
    """Archive dataset where declared P length = actual - 1 (DICOM odd-length padding).

    The length mismatch is within the ±1 tolerance and must NOT trigger a warning.
    """
    return make_hologic_archive(
        p_data=b"FAKE_P_FILE_CONTENT",
        add_padding_byte=True,
    )


@pytest.fixture
def hologic_archive_p_only_ds():
    """Archive dataset with P file data but no R file tags (partial archive)."""
    return make_hologic_archive(include_r=False)


@pytest.fixture
def hologic_archive_truncated_ds():
    """Archive where P declared length exceeds actual by 10 (triggers warning)."""
    return make_truncated_pr_archive()


@pytest.fixture
def hologic_archive_empty_ds():
    """Archive dataset with zero-length P and R blobs."""
    return make_empty_pr_archive()


@pytest.fixture
def structured_report_ds():
    """Minimal Comprehensive SR with one BMD measurement item."""
    return make_structured_report()


@pytest.fixture
def parametric_map_ds():
    """Parametric Map with FloatPixelData and RealWorldValueMappingSequence."""
    return make_parametric_map()


@pytest.fixture
def parametric_map_no_rwvm_ds():
    """Parametric Map that lacks both RWVM and float pixel data — NOT_ELIGIBLE."""
    return make_parametric_map(include_rwvm=False, include_float_pixels=False)


# ---------------------------------------------------------------------------
# File-system fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_inventory_dir(tmp_path: Path) -> Path:
    """A temporary directory pre-populated with one file of each synthetic type.

    Layout::

        <tmp_path>/
            secondary_capture.dcm
            hologic_archive.dcm
            structured_report.dcm
            parametric_map.dcm
            not_a_dicom.txt
    """
    save_to_tmpdir(make_secondary_capture(), tmp_path, "secondary_capture.dcm")
    save_to_tmpdir(make_hologic_archive(), tmp_path, "hologic_archive.dcm")
    save_to_tmpdir(make_structured_report(), tmp_path, "structured_report.dcm")
    save_to_tmpdir(make_parametric_map(), tmp_path, "parametric_map.dcm")

    txt_path = tmp_path / "not_a_dicom.txt"
    txt_path.write_text("This is not a DICOM file.\n", encoding="utf-8")

    return tmp_path
