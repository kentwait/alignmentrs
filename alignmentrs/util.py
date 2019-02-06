import os


__all__ = ['fasta_file_to_lists']


def fasta_file_to_lists(path, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str
        Location of FASTA file.
    marker_kw : str
        Keyword indicating the sample is a marker.

    Returns
    -------
    dict
        Contains list of ids, descriptions, and sequences for sample
        and marker categories.

    """
    _id = ''
    _description = ''
    _seq = ''

    sample_ids = []
    sample_descs = []
    sample_seqs = []

    marker_ids = []
    marker_descs = []
    marker_seqs = []

    if not os.path.exists(path):
        raise Exception('{} does not exist'.format(path))
    with open(path, 'r') as f:  # pylint: disable=invalid-name
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                # Store sequence if _seq has contents
                if _seq:
                    if marker_kw and (marker_kw in _id):
                        marker_ids.append(_id)
                        marker_descs.append(_description)
                        marker_seqs.append(_seq)
                    else:
                        sample_ids.append(_id)
                        sample_descs.append(_description)
                        sample_seqs.append(_seq)
                    _seq = ''
                # Split id and description
                try:
                    _id, _description = line[1:].split(' ', 1)
                except ValueError:
                    _id, _description = line[1:], ''
            else:
                _seq += line
        if _seq:
            if marker_kw and (marker_kw in _id):
                marker_ids.append(_id)
                marker_descs.append(_description)
                marker_seqs.append(_seq)
            else:
                sample_ids.append(_id)
                sample_descs.append(_description)
                sample_seqs.append(_seq)
    return {
        'sample': {
            'ids': sample_ids,
            'descriptions': sample_descs,
            'sequences': sample_seqs,
        },
        'marker': {
            'ids': marker_ids,
            'descriptions': marker_descs,
            'sequences': marker_seqs,
        }
    }
