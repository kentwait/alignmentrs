from alignmentrs.aln.classes import Alignment
from alignmentrs.aln.funcs import (fasta_file_to_alignment,
                                   mark_sites_with_chars,
                                   drop_sites_using_binary_markers)

__all__ = ['Alignment',
           'fasta_file_to_alignment', 'mark_sites_with_chars',
           'drop_sites_using_binary_markers']
