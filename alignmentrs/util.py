from collections import OrderedDict
import os
import re

from libalignmentrs.position import block_str_to_linspace
from libalignmentrs.position import simple_block_str_to_linspace


__all__ = ['fasta_file_to_lists']


def idseq_to_display(ids, chunked_sequences, template='{name}     {seq}',
                     max_length=20, id_width=15, sequence_width=55):
    if not len(ids):
        return ''
    def chunked_fn(x):
        if len(x)*len(x[0]) <= sequence_width - (len(x) - 1):
            return ' '.join(x)
        left = (sequence_width//2) // len(x[0])
        right = (sequence_width//2) // len(x[0])
        return ' '.join(x[:left]) + '...' + ' '.join(x[-right:])

    name_fn = lambda x: x + (' '*(id_width-len(x))) if len(x) <= id_width else \
        x[:id_width-3] + '...'
    seq_fn = lambda x: ''.join(x) if len(x) <= sequence_width else \
        ''.join(x[:(sequence_width//2)]) + '...' + \
        ''.join(x[-((sequence_width//2)):])

    fmt_names = (name_fn(name) for name in ids)
    fmt_seqs = ((chunked_fn(chunks) if len(chunks[0]) > 1 else seq_fn(chunks))
                if len(chunks) > 0 else ''
                for chunks in chunked_sequences)
    
    lines = [template.format(name=value[0], seq=value[1])
             for i, value in enumerate(zip(fmt_names, fmt_seqs))
             if i < max_length]
    if len(ids) > max_length:
        lines[-1] = '...'
    return '\n'.join(lines)


def add_to_history(aln, op, *args, **kwargs):
    # Add to history
    record = True
    if '_record_history' in kwargs.keys():
        record = kwargs['_record_history']
        del kwargs['_record_history']
    if record and (aln._history is not None):
        aln._history.add(op,
            args=args,
            kwargs=kwargs
        )


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


def parse_comment_list(comment_list: list):
    comments_d = dict()
    for comment in comment_list:
        k, v = comment[1:].strip().split('\t')
        if k == 'name':
            comments_d['name'] = v
        elif k == 'coords':
            comments_d['linspace'] = \
                simple_block_str_to_linspace(v)
    return comments_d


def parse_cat_comment_list(comment_list: list):
    comments_d = dict()
    subspaces = OrderedDict()
    aln_id_regex = re.compile(r'^subcoords\:(\S+)')
    for comment in comment_list:
        k, v = comment[1:].strip().split('\t')
        if k == 'name':
            comments_d['name'] = v
        elif k == 'cat_coords':
            comments_d['linspace'] = \
                block_str_to_linspace(v.lstrip('{').rstrip('}'))
        elif aln_id_regex.match(k):
            name = aln_id_regex.match(k).group(1)
            subspaces[name] = \
                simple_block_str_to_linspace(v.lstrip('{').rstrip('}'))
    comments_d['subspaces'] = subspaces
    return comments_d
