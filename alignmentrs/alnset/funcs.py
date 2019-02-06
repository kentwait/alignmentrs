from alignmentrs.alnset import AlignmentSet


__all__ = ['fasta_directory_to_alignmentset']


def fasta_directory_to_alignmentset(dirpath, name, marker_kw=None,
                                    suffix='.aln',
                                    filename_to_key_encoder=None):
    """Reads a directory containing FASTA files and stores data as a
    set of alignment objects inside an AlignmentSet.

    Parameters
    ----------
    dirpath : str
        Path containing FASTA files to be read.
    name : str
        Name of alignment set
    marker_kw : str or None, optional
        Classifies the sample as a marker if the string is
        found in the sample's ID. (default is None)
    suffix : str, optional
        Used to determine whether a file is a FASTA file (default is '.aln')
    marker_kw : str or None, optional
        Classifies the sample as a marker if the string is
        found in the sample's ID. (default is None)
    filename_to_key_encoder : function or None, optional
        If specified, the function receives the filename as input
        and outputs a key to identify a unique alignment.
        THis can be used to make sure that the same alignment
        stored as files with different filenames are not
        included multiple times.

    Returns
    -------
    AlignmentSet
        New AlignmentSet object with each FASTA file as a member Alignment
        object.

    """
    return AlignmentSet.from_fasta_dir(
        dirpath, name, marker_kw, suffix, filename_to_key_encoder)
