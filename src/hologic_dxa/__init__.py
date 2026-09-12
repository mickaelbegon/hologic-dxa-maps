"""hologic-dxa: Reproducible pipeline for Hologic DXA quantitative maps."""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("hologic-dxa")
except PackageNotFoundError:
    __version__ = "0.0.0+dev"

__all__ = ["__version__"]
