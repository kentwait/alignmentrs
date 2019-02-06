from libalignmentrs.record import Record
from libalignmentrs.alignment import BaseAlignment
from libalignmentrs import alignment as libalignment
from libalignmentrs import record as librecord

from alignmentrs.aln import Alignment, fasta_file_to_alignment
from alignmentrs.alnset import AlignmentSet, fasta_directory_to_alignmentset

from alignmentrs import aln, alnset

__author__ = 'Kent Kawashima'
__version__ = '0.6.2'
__all__ = [
    # From dynamic library
    'Record', 'BaseAlignment', 'libalignment', 'librecord',
    # Classes
    'Alignment', 'AlignmentSet',
    # Functions
    'fasta_file_to_alignment', 'fasta_directory_to_alignmentset']
