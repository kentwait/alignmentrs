from libalignmentrs import sequence as libsequence
from libalignmentrs import marker as libmarker
from libalignmentrs.sequence import Sequence
from libalignmentrs.marker import Marker
from alignmentrs.alignment import AlignmentMatrix, Alignment
from alignmentrs.alignment import fasta_file_to_list, fasta_file_to_alignment

__author__ = 'Kent Kawashima'
__version__ = '0.3.0'
__all__ = ['libalignment', 'libsequence', 'Sequence', 'Marker',
           'AlignmentMatrix', 'Alignment',
           'fasta_file_to_list', 'fasta_file_to_alignment']
