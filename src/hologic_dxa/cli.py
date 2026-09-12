"""Command-line interface for the Hologic DXA quantitative maps pipeline.

Entry point: ``hologic-dxa``  (registered in pyproject.toml project.scripts).

Each subcommand delegates to the appropriate pipeline module.  Run
``hologic-dxa --help`` or ``hologic-dxa <command> --help`` for full option
descriptions.
"""

from __future__ import annotations

import importlib.metadata
import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from hologic_dxa.logging import configure_logging, get_logger

app = typer.Typer(
    name="hologic-dxa",
    help="Reproducible pipeline for Hologic DXA quantitative maps.",
    no_args_is_help=True,
)

console = Console()
log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Shared option types
# ---------------------------------------------------------------------------

_PathArg = Annotated[Path, typer.Argument(help="Input path (file or directory).")]
_OutputOpt = Annotated[Path, typer.Option("--output", "-o", help="Output path.")]


# ---------------------------------------------------------------------------
# inspect
# ---------------------------------------------------------------------------

@app.command()
def inspect(
    path: _PathArg,
    output: _OutputOpt = Path("output/inventory"),
    log_level: Annotated[str, typer.Option("--log-level", help="Logging level.")] = "INFO",
) -> None:
    """Inspect and classify DICOM files recursively."""
    configure_logging(level=log_level)
    from hologic_dxa.dicom.inventory import build_inventory, write_inventory

    if not path.exists():
        log.error("path_not_found", path=str(path))
        console.print(f"[red]Path not found: {path}[/red]")
        raise typer.Exit(code=1)

    log.info("inspect_started", path=str(path), output=str(output))
    records = build_inventory(path)
    write_inventory(records, output)
    console.print(
        f"[green]Inventory written to {output} "
        f"({len(records)} file(s) processed).[/green]"
    )


# ---------------------------------------------------------------------------
# extract_pr
# ---------------------------------------------------------------------------

@app.command()
def extract_pr(
    path: _PathArg,
    output: _OutputOpt = Path("output/pr_files"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Extract Hologic P and R files from DICOM archives."""
    configure_logging(level=log_level)
    from hologic_dxa.dicom.extract_pr import extract_pr_from_path

    if not path.exists():
        log.error("path_not_found", path=str(path))
        console.print(f"[red]Path not found: {path}[/red]")
        raise typer.Exit(code=1)

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    extract_pr_from_path(path, output)
    console.print(f"[green]P/R files written to {output}.[/green]")


# ---------------------------------------------------------------------------
# describe_pr
# ---------------------------------------------------------------------------

@app.command()
def describe_pr(
    manifest: Annotated[Path, typer.Argument(help="Manifest file from extract_pr.")],
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Forensic description of extracted P/R files (size, entropy, magic bytes)."""
    configure_logging(level=log_level)
    from hologic_dxa.dicom.extract_pr import describe_pr

    if not manifest.exists():
        log.error("manifest_not_found", path=str(manifest))
        console.print(f"[red]File not found: {manifest}[/red]")
        raise typer.Exit(code=1)

    describe_pr(manifest)


# ---------------------------------------------------------------------------
# parse_sr
# ---------------------------------------------------------------------------

@app.command()
def parse_sr(
    path: _PathArg,
    output: _OutputOpt = Path("output/sr_results.json"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Parse Hologic APEX Structured Report DICOMs."""
    configure_logging(level=log_level)
    import json as _json
    from hologic_dxa.dicom.structured_report import load_sr_directory, sr_to_dataframe

    if not path.exists():
        log.error("path_not_found", path=str(path))
        console.print(f"[red]Path not found: {path}[/red]")
        raise typer.Exit(code=1)

    sr_records = load_sr_directory(path)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        _json.dumps([r if isinstance(r, dict) else r.__dict__ for r in sr_records], indent=2,
                    default=str),
        encoding="utf-8",
    )
    console.print(f"[green]SR results written to {output} ({len(sr_records)} record(s)).[/green]")


# ---------------------------------------------------------------------------
# normalize_sr
# ---------------------------------------------------------------------------

@app.command()
def normalize_sr(
    sr_tree: Annotated[Path, typer.Argument(help="SR results file or directory.")],
    mapping: Annotated[
        Path, typer.Option("--mapping", "-m", help="Concept mapping YAML.")
    ] = Path("config/hologic_concept_mapping.yml"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Normalize SR results using concept mapping YAML."""
    configure_logging(level=log_level)
    try:
        from hologic_dxa.dicom.sr_normalizer import normalize_sr_results
    except ImportError:
        log.error(
            "module_not_available",
            module="hologic_dxa.dicom.sr_normalizer",
            reason="not_yet_implemented",
        )
        console.print("[red]The SR normalizer module is not yet available.[/red]")
        raise typer.Exit(code=1)

    normalize_sr_results(sr_tree, mapping)


# ---------------------------------------------------------------------------
# import_maps
# ---------------------------------------------------------------------------

@app.command()
def import_maps(
    path: _PathArg,
    output: _OutputOpt = Path("output/bundle.h5"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Import pre-calibrated quantitative arrays."""
    configure_logging(level=log_level)
    from hologic_dxa.providers.exported_arrays import ExportedArraysProvider

    if not path.exists():
        log.error("path_not_found", path=str(path))
        console.print(f"[red]Path not found: {path}[/red]")
        raise typer.Exit(code=1)

    provider = ExportedArraysProvider()
    if not provider.can_read(path):
        console.print(
            f"[red]Directory '{path}' does not match the expected exported-arrays layout.\n"
            "Required files: fat.npy (or fat.csv), lean.npy, bmc.npy, metadata.json.[/red]"
        )
        raise typer.Exit(code=1)

    bundle = provider.load(path)
    output = Path(output)
    _write_simple_bundle_hdf5(bundle, output)
    console.print(f"[green]Bundle written to {output}.[/green]")


# ---------------------------------------------------------------------------
# validate_maps
# ---------------------------------------------------------------------------

@app.command()
def validate_maps(
    bundle: Annotated[Path, typer.Argument(help="Path to bundle HDF5 file.")],
    sr: Annotated[
        Path | None, typer.Option("--sr", help="SR results JSON or CSV.")
    ] = None,
    apex_csv: Annotated[
        Path | None, typer.Option("--apex-csv", help="APEX CSV export.")
    ] = None,
    output: _OutputOpt = Path("output/validation"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Validate mass conservation against SR/APEX reference values."""
    configure_logging(level=log_level)
    import pandas as pd
    from hologic_dxa.validation.mass_conservation import validate_mass_conservation

    bundle_obj = _load_bundle(bundle)

    sr_df: pd.DataFrame | None = None
    if sr is not None:
        sr_df = _load_dataframe(sr, "SR results")

    apex_df: pd.DataFrame | None = None
    if apex_csv is not None:
        apex_df = _load_dataframe(apex_csv, "APEX CSV")

    report = validate_mass_conservation(bundle_obj, sr_results=sr_df, apex_csv=apex_df)

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "validation_report.json"
    json_path.write_text(report.to_json(), encoding="utf-8")
    csv_path = output / "validation_report.csv"
    report.to_dataframe().to_csv(csv_path, index=False)

    status = "[green]PASSED[/green]" if report.overall_passed else "[red]FAILED[/red]"
    console.print(f"Validation: {status}  ({len(report.comparisons)} comparison(s))")
    for w in report.warnings:
        console.print(f"  [yellow]WARNING:[/yellow] {w}")
    console.print(f"Report written to {output}.")
    if not report.overall_passed:
        raise typer.Exit(code=1)


# ---------------------------------------------------------------------------
# compute_mass_properties
# ---------------------------------------------------------------------------

@app.command()
def compute_mass_properties(
    bundle: Annotated[Path, typer.Argument(help="Path to bundle HDF5 file.")],
    output: _OutputOpt = Path("output/mass_properties.json"),
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Compute mass, CoM, and inertia from quantitative maps."""
    configure_logging(level=log_level)
    import json

    from hologic_dxa.maps.integration import compute_mass_properties as _compute

    bundle_obj = _load_bundle(bundle)
    props = _compute(bundle_obj)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(props, indent=2), encoding="utf-8")
    console.print(
        f"[green]Mass properties written to {output} "
        f"({len(props)} region(s)).[/green]"
    )


# ---------------------------------------------------------------------------
# export_hdf5
# ---------------------------------------------------------------------------

@app.command()
def export_hdf5(
    bundle: Annotated[Path, typer.Argument(help="Path to bundle HDF5 file.")],
    output: _OutputOpt = Path("output/bundle_export.h5"),
    compression: Annotated[
        str, typer.Option("--compression", help="HDF5 compression filter.")
    ] = "gzip",
    compression_level: Annotated[
        int, typer.Option("--compression-level", help="Compression level (0–9).")
    ] = 4,
    log_level: Annotated[str, typer.Option("--log-level")] = "INFO",
) -> None:
    """Export map bundle to HDF5."""
    configure_logging(level=log_level)
    from hologic_dxa.export.hdf5 import export_bundle_hdf5

    bundle_obj = _load_bundle(bundle)
    export_bundle_hdf5(
        bundle_obj,
        Path(output),
        compression=compression,
        compression_opts=compression_level,
    )
    console.print(f"[green]Bundle exported to {output}.[/green]")


# ---------------------------------------------------------------------------
# doctor
# ---------------------------------------------------------------------------

@app.command()
def doctor() -> None:
    """Report pipeline status: installed deps, available providers, blocked features."""
    _run_doctor()


def _run_doctor() -> None:
    """Implementation of the doctor command, separated for testability."""
    _check_python()
    _check_packages()
    _check_providers()
    _check_features()


def _check_python() -> None:
    table = Table(title="Python Runtime", show_header=True, header_style="bold cyan")
    table.add_column("Item", style="bold")
    table.add_column("Value")
    table.add_column("Status")

    py_ver = sys.version.split()[0]
    py_ok = tuple(int(x) for x in py_ver.split(".")[:2]) >= (3, 11)
    table.add_row(
        "Python version",
        py_ver,
        "[green]OK[/green]" if py_ok else "[red]REQUIRES >=3.11[/red]",
    )

    pkg_ver = _safe_version("hologic-dxa")
    table.add_row("hologic-dxa", pkg_ver, "[green]OK[/green]" if pkg_ver != "?" else "[yellow]?[/yellow]")

    console.print(table)


def _check_packages() -> None:
    table = Table(title="Key Dependencies", show_header=True, header_style="bold cyan")
    table.add_column("Package", style="bold")
    table.add_column("Version")
    table.add_column("Status")

    packages = [
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("pandas", "pandas"),
        ("pydicom", "pydicom"),
        ("highdicom", "highdicom"),
        ("pydantic", "pydantic"),
        ("h5py", "h5py"),
        ("pint", "pint"),
        ("structlog", "structlog"),
        ("rich", "rich"),
        ("typer", "typer"),
        ("pyyaml", "PyYAML"),
    ]

    for import_name, pkg_name in packages:
        ver = _safe_version(pkg_name)
        try:
            __import__(import_name)
            status = "[green]OK[/green]"
        except ImportError:
            status = "[red]MISSING[/red]"
            ver = "not installed"
        table.add_row(pkg_name, ver, status)

    console.print(table)


def _check_providers() -> None:
    table = Table(title="Quantitative Map Providers", show_header=True, header_style="bold cyan")
    table.add_column("Provider", style="bold")
    table.add_column("Status")
    table.add_column("Notes")

    providers = [
        (
            "Pre-calibrated array import",
            _module_available("hologic_dxa.maps.importer"),
            "Import quantitative arrays exported from third-party software.",
        ),
        (
            "Hologic SDK",
            False,
            "Requires proprietary Hologic SDK (not open-source). "
            "Set HOLOGIC_SDK_PATH to enable.",
        ),
        (
            "DICOM Parametric Map",
            True,
            "Standard DICOM PM objects are supported via pydicom/highdicom.",
        ),
        (
            "SR result extraction",
            _module_available("hologic_dxa.dicom.sr_parser"),
            "Parse regional results from Hologic APEX Structured Reports.",
        ),
    ]

    for name, available, notes in providers:
        status = "[green]Available[/green]" if available else "[red]Unavailable[/red]"
        table.add_row(name, status, notes)

    console.print(table)


def _check_features() -> None:
    table = Table(
        title="Pipeline Features",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Feature", style="bold")
    table.add_column("Status")
    table.add_column("Reason / Notes")

    operational = [
        ("DICOM inventory & classification", True, ""),
        ("HDF5 bundle export", True, ""),
        ("Mass / CoM / inertia computation", True, ""),
        ("Mass conservation validation", True, "Requires reference SR or APEX CSV."),
        ("Unit conversion (g/cm2 ↔ mg/cm2)", True, ""),
        ("QC: total consistency check", True, ""),
    ]

    blocked = [
        (
            "P/R file extraction",
            False,
            "Requires pr_extractor module (not yet implemented).",
        ),
        (
            "SR parsing",
            False,
            "Requires sr_parser module (not yet implemented).",
        ),
        (
            "Map import from array files",
            False,
            "Requires maps.importer module (not yet implemented).",
        ),
        (
            "BodyLoop integration",
            False,
            "Requires trimesh/open3d and bodyloop module.",
        ),
        (
            "Hologic SDK calibration",
            False,
            "Requires proprietary Hologic SDK (not distributable).",
        ),
    ]

    for name, status, notes in operational + blocked:
        symbol = "[green]Operational[/green]" if status else "[red]Blocked[/red]"
        table.add_row(name, symbol, notes)

    console.print(table)


# ---------------------------------------------------------------------------
# Internal helpers shared across commands
# ---------------------------------------------------------------------------

def _safe_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "?"


def _module_available(module: str) -> bool:
    try:
        __import__(module)
        return True
    except ImportError:
        return False


def _load_bundle(path: Path) -> object:
    """Load a QuantitativeMapBundle from an HDF5 file.

    Returns a lightweight namespace populated with numpy arrays read from the
    standard HDF5 structure written by :func:`_write_simple_bundle_hdf5` or
    :func:`export_bundle_hdf5`.
    """
    import h5py
    import types

    path = Path(path)
    if not path.exists():
        log.error("bundle_not_found", path=str(path))
        console.print(f"[red]Bundle file not found: {path}[/red]")
        raise typer.Exit(code=1)

    try:
        with h5py.File(path, "r") as f:
            bundle = types.SimpleNamespace()

            # Maps: fat, lean, bmc, total — each exposed as .values
            maps_ns = types.SimpleNamespace()
            for component in ("fat", "lean", "bmc", "total"):
                ds_name = f"maps/{component}_areal_density_g_cm2"
                arr = f[ds_name][:] if ds_name in f else None
                qmap_ns = types.SimpleNamespace(
                    values=arr,
                    pixel_area_cm2=None,
                    x_coordinates_mm=None,
                    y_coordinates_mm=None,
                )
                setattr(maps_ns, f"{component}_areal_density_g_cm2", arr)
                setattr(bundle, component, qmap_ns if arr is not None else None)
            bundle.maps = maps_ns

            # Masks
            masks_ns = types.SimpleNamespace(
                body=f["masks/body"][:].astype(bool) if "masks/body" in f else None,
                bone=f["masks/bone"][:].astype(bool) if "masks/bone" in f else None,
                valid=f["masks/valid"][:].astype(bool) if "masks/valid" in f else None,
            )
            bundle.masks = masks_ns
            bundle.body_mask = masks_ns.body
            bundle.bone_mask = masks_ns.bone
            bundle.valid_mask = (
                masks_ns.valid if masks_ns.valid is not None else masks_ns.body
            )
            if bundle.valid_mask is None:
                # Fall back: all pixels valid
                first_arr = next(
                    (
                        getattr(bundle, c).values
                        for c in ("fat", "lean", "bmc", "total")
                        if getattr(bundle, c) is not None
                    ),
                    None,
                )
                if first_arr is not None:
                    bundle.valid_mask = np.ones(first_arr.shape, dtype=bool)

            # Geometry
            pixel_area = (
                float(f["geometry/pixel_area_cm2"][...])
                if "geometry/pixel_area_cm2" in f
                else None
            )
            x_mm = f["geometry/x_coordinates_mm"][:] if "geometry/x_coordinates_mm" in f else None
            y_mm = f["geometry/y_coordinates_mm"][:] if "geometry/y_coordinates_mm" in f else None

            if x_mm is not None:
                for c in ("fat", "lean", "bmc", "total"):
                    qm = getattr(bundle, c)
                    if qm is not None:
                        qm.x_coordinates_mm = x_mm
                        qm.y_coordinates_mm = y_mm
                        qm.pixel_area_cm2 = pixel_area

            geom_ns = types.SimpleNamespace()
            geom_ns.pixel_area_cm2 = pixel_area

            def _pixel_area_cm2_scalar() -> float:
                if pixel_area is None:
                    raise ValueError(
                        "pixel_area_cm2 is missing from the HDF5 bundle geometry group."
                    )
                return pixel_area

            geom_ns.pixel_area_cm2_scalar = _pixel_area_cm2_scalar  # type: ignore[attr-defined]
            # Reconstruct geometry shape from arrays
            if x_mm is not None:
                geom_ns.rows = x_mm.shape[0]
                geom_ns.columns = x_mm.shape[1]
                geom_ns.pixel_spacing_mm = (1.0, 1.0)  # unknown; coordinates already stored
                geom_ns.origin_mm = (float(x_mm[0, 0]), float(y_mm[0, 0]))  # type: ignore[index]
            bundle.geometry = geom_ns
            bundle.x_coordinates_mm = x_mm
            bundle.y_coordinates_mm = y_mm

            # Metadata / provenance
            meta_attrs = dict(f["metadata"].attrs) if "metadata" in f else {}
            bundle.metadata = meta_attrs
            bundle.provenance = None

    except OSError as exc:
        log.error("bundle_read_error", path=str(path), error=str(exc))
        console.print(f"[red]Cannot read bundle file '{path}': {exc}[/red]")
        raise typer.Exit(code=1) from exc

    return bundle


def _write_simple_bundle_hdf5(bundle: object, output_path: Path) -> None:
    """Write a ``maps.bundle.QuantitativeMapBundle`` to HDF5.

    Produces the same structure as :func:`export_bundle_hdf5` so that
    downstream commands can read it back via :func:`_load_bundle`.
    """
    import h5py
    import numpy as _np

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fat: _np.ndarray = getattr(bundle, "fat")
    lean: _np.ndarray = getattr(bundle, "lean")
    bmc: _np.ndarray = getattr(bundle, "bmc")
    total: _np.ndarray = fat + lean + bmc
    valid_mask: _np.ndarray = getattr(bundle, "valid_mask")
    pixel_area: float = getattr(bundle, "pixel_area_cm2")

    with h5py.File(output_path, "w") as f:
        maps_grp = f.create_group("maps")
        for name, arr in (
            ("fat_areal_density_g_cm2", fat),
            ("lean_areal_density_g_cm2", lean),
            ("bmc_areal_density_g_cm2", bmc),
            ("total_areal_density_g_cm2", total),
        ):
            ds = maps_grp.create_dataset(
                name, data=arr.astype(_np.float32), compression="gzip", compression_opts=4
            )
            ds.attrs["units"] = "g/cm^2"

        masks_grp = f.create_group("masks")
        ds_valid = masks_grp.create_dataset(
            "valid", data=valid_mask.astype(_np.uint8), compression="gzip", compression_opts=4
        )
        ds_valid.attrs["units"] = "boolean (0=False, 1=True)"

        geom_grp = f.create_group("geometry")
        ds_area = geom_grp.create_dataset("pixel_area_cm2", data=_np.float64(pixel_area))
        ds_area.attrs["units"] = "cm^2"

        meta_grp = f.create_group("metadata")
        meta_grp.attrs["source_sop_instance_uid"] = str(
            getattr(bundle, "source_sop_instance_uid", "")
        )
        meta_grp.attrs["provider_name"] = str(getattr(bundle, "provider_name", ""))
        meta_grp.attrs["already_geometry_corrected"] = bool(
            getattr(bundle, "already_geometry_corrected", False)
        )

        f.create_group("provenance")

    log.info("simple_bundle_hdf5_written", path=str(output_path))


def _load_dataframe(path: Path, label: str) -> object:
    """Load a CSV or JSON file as a pandas DataFrame."""
    import pandas as pd

    path = Path(path)
    if not path.exists():
        log.error("file_not_found", path=str(path), label=label)
        console.print(f"[red]{label} file not found: {path}[/red]")
        raise typer.Exit(code=1)

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in (".json", ".jsonl"):
        return pd.read_json(path)
    log.error("unsupported_file_format", path=str(path), suffix=suffix)
    console.print(
        f"[red]Unsupported format '{suffix}' for {label}. Use .csv or .json.[/red]"
    )
    raise typer.Exit(code=1)
