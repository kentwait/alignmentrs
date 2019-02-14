from libalignmentrs.alignment import fasta_file_to_basealignments
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
    # Create alignments
    if marker_kw is None:
        marker_kw = ''
    return Alignment(name, *fasta_file_to_basealignments(path, marker_kw))


# def split_concatenated_alignment(aln, catblocks=None,
#                                  description_decoder=None):
#     pass

# def blocks_list_to_df(blocks_list):
#     pass

# def catblocks_list_to_df(catblocks_list):
#     pass
