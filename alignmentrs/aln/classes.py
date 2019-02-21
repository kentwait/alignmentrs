from collections import OrderedDict
import itertools
import os
from copy import deepcopy

import pandas as pd

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import BaseRecord
from libalignmentrs.readers import fasta_to_records


__all__ = ['Alignment', 'CatAlignment']


class _Rows:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 0

    def insert(self, position, records, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # TODO: Check data type of position
        if isinstance(records, BaseRecord):
            aln._alignment.insert_rows(position, records)
        elif isinstance(records, list) and \
            sum((isinstance(rec, BaseAlignment) for rec in records)):
            aln._alignment.insert_rows(position, records)
        else:
            raise TypeError('records must be a BaseRecord or a list of BaseRecord objects')
        if copy is True:
            return aln

    def prepend(self, records, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(records, BaseRecord):
            aln._alignment.insert_rows(0, records)
        elif isinstance(records, list) and \
            sum((isinstance(rec, BaseAlignment) for rec in records)):
            aln._alignment.insert_row(0, records)
        else:
            raise TypeError('records must be a BaseRecord or a list of BaseRecord objects')
        if copy is True:
            return aln

    def append(self, records, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(records, BaseRecord):
            aln._alignment.append_rows(records)
        elif isinstance(records, list) and \
            sum((isinstance(rec, BaseAlignment) for rec in records)):
            aln._alignment.append_row(records)
        else:
            raise TypeError('records must be a BaseRecord or a list of BaseRecord objects')
        if copy is True:
            return aln

    def remove(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            aln._alignment.remove_row(positions)
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            aln._alignment.remove_rows(positions)
        else:
            raise TypeError('positions must be an int or a list of int')
        if copy is True:
            return aln

    def retain(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            aln._alignment.retain_row(positions)
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            aln._alignment.retain_rows(positions)
        else:
            raise TypeError('positions must be an int or a list of int')
        if copy is True:
            return aln

    def drain(self, positions):
        raise NotImplementedError()

    def replace(self, positions, records, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # TODO: Check data type of positions
        if isinstance(records, BaseRecord):
            aln._alignment.replace_row(positions, records)
        elif isinstance(records, list) and \
            sum((isinstance(rec, BaseAlignment) for rec in records)):
            aln._alignment.replace_rows(positions, records)
        else:
            raise TypeError('records must be a BaseRecord or a list of BaseRecord objects')
        if copy is True:
            return aln

    def reorder(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            aln._alignment.reorder_rows([positions])
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            aln._alignment.reorder_rows(positions)
        else:        
            raise TypeError('positions must be an int or a list of int')
        if copy is True:
            return aln

    def filter(self, function, copy=False, dry_run=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # Function accepts a Record, outputs true or false
        if not(function is not None and callable(function)):
            raise TypeError('missing filter function')
        positions = [i for i in range(aln.nrows) 
                     if function(aln._alignment.get_record(i))]
        remove_positions = aln._alignment.invert_rows(positions)
        if dry_run:
            parts = []
            parts.append('[Filter]')
            parts.append('True = {}/{}'.format(
                len(positions), aln.nrows))
            parts.append('False = {}/{}'.format(
                len(remove_positions), aln.nrows))
            print('\n'.join(parts))
            return {'function': function, True: positions, False: remove_positions}
        aln.rows.remove(remove_positions)
        if copy is True:
            return aln

    def iter(self):
        for i in range(self._instance.nrows):
            yield self._instance._alignment.get_record(i)

    def iter_sequences(self):
        for i in range(self._instance.nrows):
            yield self._instance._alignment.get_row(i)

    def __iter__(self):
        return self.iter()

    def __len__(self):
        return self._instance.nrows


class _Cols:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 1

    def insert(self, position, values, copy=False):
        raise NotImplementedError()

    def prepend(self, values, copy=False):
        raise NotImplementedError()

    def append(self, values, copy=False):
        raise NotImplementedError()

    def remove(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')
        retain_positions = aln._alignment.invert_cols(positions)        
        aln._alignment.remove_cols(positions)
        aln._column_metadata = aln._column_metadata.iloc[retain_positions]
        if copy is True:
            return aln

    def retain(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            pass
        else:        
            raise TypeError('positions must be an int or a list of int')
        remove_positions = aln._alignment.invert_cols(positions)
        aln._alignment.remove_cols(remove_positions)
        aln._column_metadata = aln._column_metadata.iloc[positions]
        if copy is True:
            return aln

    def drain(self, positions):
        raise NotImplementedError()

    def replace(self, positions, values, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # TODO: Check data type of positions
        if isinstance(values, list):
            if sum((isinstance(val, str) for val in values)):
                values = [[val] for val in values]
            elif sum((isinstance(val, list) for val in values)) and \
                sum((isinstance(val, str) for lst in values for val in lst)):
                pass
            else:
                raise TypeError('values must be a list of str, or list of list of str')
        else:
            raise TypeError('values must be a list of str, or list of list of str')
        aln._alignment.replace_cols(positions, values)
        if copy is True:
            return aln

    def reorder(self, positions, copy=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')
        aln._alignment.reorder_cols(positions)
        aln._column_metadata = aln._column_metadata.iloc[positions]
        if copy is True:
            return aln

    def filter(self, function, copy=False, dry_run=False):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # Function accepts a list of str, outputs true or false
        if not(function is not None and callable(function)):
            raise TypeError('missing filter function')
        positions = [i for i in range(aln.ncols) 
                     if function(aln._alignment.get_col(i))]
        remove_positions = aln._alignment.invert_cols(positions)
        if dry_run:
            parts = []
            parts.append('[Filter]')
            parts.append('True = {}/{}'.format(
                len(positions), aln.ncols))
            parts.append('False = {}/{}'.format(
                len(remove_positions), aln.ncols))
            print('\n'.join(parts))
            return {'function': function, True: positions, False: remove_positions}
        aln.cols.remove(remove_positions)
        if copy is True:
            return aln

    def iter(self, skip_n=None, chunk_size=None):
        cnt = 0
        out = []
        if skip_n is None:
            if chunk_size is None:
                skip_n = 1
            else:
                skip_n = chunk_size
        if chunk_size is None:
            chunk_size = 1

        for i in range(0, self._instance.ncols-(chunk_size-1), skip_n):
            if chunk_size == 1:
                yield self._instance._alignment.get_col(i)
            else:
                yield self._instance._alignment.get_chunk(i, chunk_size)

    def __iter__(self):
        return self.iter()
    
    def __len__(self):
        return self._instance.ncols


class Alignment:
    """Reperesents a multiple sequence alignment of samples.

    The Alignment object encapsulates information generally
    included in the FASTA format:
    - sequence names/ids
    - descriptions
    - sequences

    Additionally, the Alignment object also stores comments
    (lines prefixed by a semicolon ";") as metadata.

    Attributes
    ----------
    name : str
        Name of the alignment.
    samples : BaseAlignment
        Alignment of sample sequences.
    metadata : dict
        Other information related to the alignment.
    """
    def __init__(self, name, records, chunk_size: int=1,
                 index=None, metadata: dict=None, column_metadata=None,
                 **kwargs):
                # TODO: Update docstrings
        """Creates a new Alignment object from a sample BaseAlignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        sample_alignment : BaseAlignment
            Alignment of sample sequences.
        linspace : BlockSpace, optional
            Linear space that assigns coordinate values to alignment
            columns. (default is None, the constructor will create a new
            linear space starting from 0).
        metadata : dict, optional
            Other information related to the alignment. (default is None,
            which creates a blank dictionary)
        **kwargs
            Other keyword arguments used to initialize states in the
            Alignment object.

        Raises
        ------
        ValueError
            Alignment is instantiated with an empty sample alignment, or
            instantiated with sample and marker alingments of unequal number of
            sites.

        """
        self.name = name
        self._chunk_size = chunk_size
        self._alignment: BaseAlignment = \
            self._alignment_constructor(records, chunk_size)
        self._index = self._index_constructor(index)
        self.metadata = self._metadata_constructor(metadata)
        self._column_metadata = \
            self._col_metadata_constructor(column_metadata, self.index)
        self._rows = _Rows(self)
        self._cols = _Cols(self)

    # Constructors
    def _alignment_constructor(self, records, chunk_size):
        if isinstance(records, list):
            if not sum((isinstance(rec, BaseRecord) for rec in records)):
                raise TypeError('records must be a list of BaseRecord objects')
            return BaseAlignment(records, chunk_size)
        elif isinstance(records, BaseAlignment):
            if records.chunk_size != chunk_size:
                records.chunk_size = chunk_size
            return records
        raise TypeError('records must be a list of BaseRecord objects or a BaseAlignment')

    def _metadata_constructor(self, metadata):
        if metadata is None:
            return dict()
        elif isinstance(metadata, dict):
            return metadata
        raise TypeError('metadata must be a dictionary object')

    def _index_constructor(self, index):
        if index is None:
            return pd.Index(list(range(self._alignment.ncols)))
        elif isinstance(index, pd.Index):
            return index
        raise TypeError(
            'index must be a {} object'.format(pd.Index.__mro__[0]))

    def _col_metadata_constructor(self, column_metadata, index):
        if column_metadata is None:
            return pd.DataFrame(None, index=index)
        elif isinstance(column_metadata, dict):
            # Check if values match the length of the index
            for key, val in column_metadata.items():
                if len(val) != len(self.index):
                    raise ValueError('{} value length does not match the number of columns'.format(key))
            return pd.DataFrame(column_metadata, index=self.index)
        elif isinstance(column_metadata, pd.DataFrame):
            if len(column_metadata) != len(self.index):
                raise ValueError('length of column_metadata dataframe does not match the number of columns'.format(key))
            return column_metadata
        raise TypeError('column_metadata must be a dictionary or a {} object'.format(pd.DataFrame.__mro__[0]))

    # Properties
    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    @property
    def index(self):
        """pandas.core.indexes.base.Index: Returns the column index
        of the alignment."""
        return self._index

    @property
    def chunk_size(self):
        """int: Returns the chunk size of the alignment."""
        return self._chunk_size

    @property
    def column_metadata(self):
        """pandas.core.frame.DataFrame: Returns the associated column
        metadata as a pandas DataFrame."""
        return self._column_metadata

    @property
    def nrows(self):
        """int: Returns the number of rows in the alignment."""
        return self._alignment.nrows

    @property
    def ncols(self):
        """int: Returns the number of columns in the alignment."""
        return self._alignment.ncols

    @property
    def nchars(self):
        """int: Returns the number of aligned characters
        (ncols * chunk_size)."""
        return self._alignment.nchars

    @property
    def ids(self):
        """list of str: Returns the list of identifiers."""
        return self._alignment.ids

    @property
    def descriptions(self):
        """list of str: Returns the list of descriptions."""
        return self._alignment.descriptions

    @property
    def sequences(self):
        """list of str: Returns the list of sequences."""
        return self._alignment.sequences

    @property
    def chunked_sequences(self):
        """list of list of str: Returns the list of sequences."""
        return self._alignment.chunked_sequences

    # Methods
    # ==========================================================================

    def to_records(self):
        return self._alignment.get_records(list(range(self.nrows)))

    def copy(self):
        return self.__class__(
            self.name, self._alignment.copy(),
            chunk_size=self.chunk_size,
            index=self._index.copy(deep=True),
            metadata=deepcopy(self.metadata),
            column_metadata=self._column_metadata.copy(deep=True)
        )

    def set_chunk_size(self, value, func_map=None, copy=False):
        aln = self
        if copy is True:
            aln = self.copy()
        aln._alignment.chunk_size = value
        aln._chunk_size = value
        if copy is True:
            return aln
        # TODO: add a way for the positional data to adjust

    # TODO: implement __copy__ and __deepcopy__

    # Special methods
    # ==========================================================================

    def __getitem__(self, key):
        if isinstance(key, str):
            if key in self.ids:
                i = self._alignment.row_names_to_indices([key])
                return self.samples.get_row(i[0])
            raise KeyError('key did not match any identifier')
        elif isinstance(key, int):
            return self.samples.get_col(key)
        raise TypeError('key must be str or int')

    def __delitem__(self, key):
        if isinstance(key, str):
            if key in self.ids:
                i = self._alignment.row_names_to_indices([key])
                return self._alignment.remove_rows(i)
            raise KeyError('key did not match any identifier')
        elif isinstance(key, int):
            return self.remove_cols(key)
        raise TypeError('key must be str or int')

    def __iter__(self):
        raise NotImplementedError(
            'Use .rows.iter() or .rows.iter_sequences() '
            'to iterate over records or sample sequences, respectively.\n'
            'Use .cols.iter() to iterate across columns in the alignment.')

    def __repr__(self):
        parts = []
        parts.append('[Alignment]')
        parts.append('name = {}'.format(self.name))
        parts.append('nrows = {}'.format(self.nrows))
        parts.append('ncols = {}'.format(self.ncols))
        parts.append('chunk_size = {}'.format(self.chunk_size))
        if self:
            aln = idseq_to_display(self.ids, self.chunked_sequences)
            parts += ['', aln, '']
        parts.append('metadata_keys = [{}]'.format(
            ', '.join(list(self.metadata.keys()))
        ))
        parts.append('column_metadata_keys = [{}]'.format(
            ', '.join(list(self._column_metadata.keys()))
        ))
        return '\n'.join(parts)

    def __str__(self):
        return str(self._alignment)

    def __len__(self):
        raise NotImplementedError(
            'Use .nrows to get the number of samples, or '
            '.ncols to get the number of columns in the alignment.')

    def __bool__(self):
        if self.ncols == 0 or self.nrows == 0:
            return False
        return True

    def __hash__(self):
        return hash((
            str(self),
            self._linspace.to_block_str(),
            tuple(self.metadata.items())
        ))

    def __eq__(self, other):
        return hash(self) == hash(other)


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
