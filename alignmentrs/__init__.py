from libalignmentrs.alignment import BaseAlignment, fasta_file_to_basealignments
from libalignmentrs.record import Record, fasta_file_to_records
from libalignmentrs.position import Block, LSpace
import libalignmentrs as librs

from alignmentrs import aln, alnset
from alignmentrs.aln import Alignment, fasta_file_to_alignment
from alignmentrs.alnset import AlignmentSet, fasta_directory_to_alignmentset


__author__ = 'Kent Kawashima'
__version__ = '0.8.3'
__all__ = [
    # From dynamic library
    'BaseAlignment', 'Record', 'Block', 'LSpace', 'librs',
    # Modules
    'aln', 'alnset',
    # Classes
    'Alignment', 'AlignmentSet',
    # Functions
    'fasta_file_to_alignment', 'fasta_directory_to_alignmentset',
    'fasta_file_to_records']
