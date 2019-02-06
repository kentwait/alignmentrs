from libalignmentrs.alignment import BaseAlignment
from alignmentrs.util import fasta_file_to_lists
from alignmentrs.aln import Alignment


__all__ = ['fasta_file_to_alignment']


def fasta_file_to_alignment(path, name, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str
        Location of FASTA file.
    name : str
        Name of the alignment.
    marker_kw : str
        Keyword indicating the sample is a marker.

    Returns
    -------
    Alignment

    """
    d = fasta_file_to_lists(path, marker_kw=marker_kw)
    sample_aln = BaseAlignment(d['sample']['ids'],
                               d['sample']['descriptions'],
                               d['sample']['sequences'])
    if len(d['marker']['ids']) > 0:
        marker_aln = BaseAlignment(d['marker']['ids'],
                                   d['marker']['descriptions'],
                                   d['marker']['sequences'])
    else:
        marker_aln = None
    # Create alignments
    return Alignment(name, sample_aln, marker_aln)

def split_concatenated_alignment(aln, catblocks=None,
                                 description_decoder=None):
    pass

def blocks_list_to_df(blocks_list):
    pass

def catblocks_list_to_df(catblocks_list):
    pass
