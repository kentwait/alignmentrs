from libalignmentrs.record import Record
from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.fasta import fasta_file_to_basealignments
import libalignmentrs as librs

from alignmentrs import aln, alnset
from alignmentrs.aln import Alignment, fasta_file_to_alignment
from alignmentrs.alnset import AlignmentSet, fasta_directory_to_alignmentset


__author__ = 'Kent Kawashima'
__version__ = '0.7.2'
__all__ = [
    # From dynamic library
    'Record', 'BaseAlignment', 'librs',
    # Modules
    'aln', 'alnset',
    # Classes
    'Alignment', 'AlignmentSet',
    # Functions
    'fasta_file_to_alignment', 'fasta_directory_to_alignmentset']
