import pkg_resources

__version__ = pkg_resources.get_distribution("bamread").version

from bamread.read import read_bam, read_bam_full

read_sam = read_bam
read_sam_full = read_bam_full
