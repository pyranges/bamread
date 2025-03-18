try:
    # For Python 3.8+
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # For older Python versions, install the backport:
    # pip install importlib_metadata
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("bamread")
except PackageNotFoundError:
    __version__ = "unknown"

from bamread.read import read_bam, read_bam_full

read_sam = read_bam
read_sam_full = read_bam_full
