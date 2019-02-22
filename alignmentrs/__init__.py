from libalignmentrs.readers import fasta_to_records
import libalignmentrs as librs

from alignmentrs import aln
# from alignmentrs import alnset
from alignmentrs.aln import Alignment
# from alignmentrs.alnset import AlignmentSet, fasta_directory_to_alignmentset


__author__ = 'Kent Kawashima'
__version__ = '0.9.0'
__all__ = [
    # From dynamic library
    'librs',
    # Modules
    'aln', 'alnset',
    # Classes
    'Alignment', 'AlignmentSet',
    # Functions
    'fasta_file_to_alignment', 'fasta_directory_to_alignmentset',
    'fasta_to_records']
