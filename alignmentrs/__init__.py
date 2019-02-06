from libalignmentrs.record import Record
from libalignmentrs.alignment import BaseAlignment
from alignmentrs.alignment import Alignment
from alignmentrs.alignment import fasta_file_to_alignment
from alignmentrs.set import AlignmentSet
import libalignmentrs.record as librecord


__author__ = 'Kent Kawashima'
__version__ = '0.6.2'
__all__ = ['librecord', 'Record', 'BaseAlignment',
           'Alignment', 'AlignmentSet',
           'fasta_file_to_alignment']
