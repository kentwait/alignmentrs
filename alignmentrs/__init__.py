from libalignmentrs.record import Record
from libalignmentrs.alignment import BaseAlignment
from libalignmentrs import alignment as libalignment
from libalignmentrs import record as librecord

from alignmentrs.alignment import Alignment, fasta_file_to_alignment
from alignmentrs.set import AlignmentSet
from alignmentrs import alignment, alnset

__author__ = 'Kent Kawashima'
__version__ = '0.6.2'
__all__ = ['Record', 'BaseAlignment','Alignment', 'AlignmentSet',
           'libalignment', 'librecord', 'alignment', 'alnset',
           'fasta_file_to_alignment']
