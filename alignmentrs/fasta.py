from alignmentrs import Sequence, Marker
from alignmentrs.alignment import Alignment

def fasta_file_to_list(path, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str

    Returns
    -------
    list of tuple

    """
    name = ''
    description = ''
    _seq = ''
    seq_list = []
    with open(path, 'r') as f:  # pylint: disable=invalid-name
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                # Store sequence if _seq has contents
                if _seq:
                    if marker_kw:
                        if marker_kw in name:
                            seq = Marker(name, description, _seq)
                        else:
                            seq = Sequence(name, description, _seq)
                    else:
                        seq = Sequence(name, description, _seq)
                    seq_list.append(seq)
                    _seq = ''
                # Split id and description
                try:
                    name, description = line[1:].split(' ', 1)
                except ValueError:
                    name = line[1:]
                    description = ''
            else:
                _seq += line
        if _seq:
            if marker_kw is not None:
                if marker_kw in name:
                    seq = Marker(name, description, _seq)
                else:
                    seq = Sequence(name, description, _seq)
            else:
                seq = Sequence(name, description, _seq)
            seq_list.append(seq)
    return seq_list

def fasta_file_to_alignment(path, marker_kw=None,
                            sample_to_uint_fn=None, uint_to_sample_fn=None,
                            marker_to_uint_fn=None, uint_to_marker_fn=None):
    """Reads a FASTA formatted text file into an Alignment object.

    Parameters
    ----------
    path : str

    Returns
    -------
    list of tuple

    """
    return Alignment.from_fasta(path, marker_kw,
                                sample_to_uint_fn, uint_to_sample_fn,
                                marker_to_uint_fn, uint_to_marker_fn)
