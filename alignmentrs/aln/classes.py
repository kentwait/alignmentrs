from collections import Counter
import os
from copy import deepcopy

import pandas
import numpy

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import BaseRecord
from alignmentrs.util import idseq_to_display
from alignmentrs.aln.mixins import (RecordsSerdeMixin, FastaSerdeMixin,
                                    JsonSerdeMixin, PickleSerdeMixin)
from .mutator import RowMutator, ColMutator


__all__ = ['Alignment', 'CatAlignment']


class Alignment(PickleSerdeMixin, JsonSerdeMixin, FastaSerdeMixin, 
                RecordsSerdeMixin, object):
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
        self._alignment: BaseAlignment = \
            self._alignment_constructor(records, chunk_size)
        self._index = self._index_constructor(index)
        self.metadata = self._metadata_constructor(metadata)
        self._column_metadata = \
            self._col_metadata_constructor(column_metadata, self.index)
        self._rows = RowMutator(self)
        self._cols = ColMutator(self)
        # TODO: Add logging to log and replay actions performed using the API

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
            return pandas.Index(range(self._alignment.ncols))
        elif isinstance(index, pandas.Index):
            return index
        raise TypeError(
            'index must be a {} object'.format(pandas.Index.__mro__[0]))

    def _col_metadata_constructor(self, column_metadata, index):
        if column_metadata is None:
            df = pandas.DataFrame(None, index=index)
        elif isinstance(column_metadata, dict):
            # Check if values match the length of the index
            for key, val in column_metadata.items():
                if len(val) != len(self.index):
                    raise ValueError('{} value length does not match the number of columns'.format(key))
            df = pandas.DataFrame(column_metadata, index=self.index)
        elif isinstance(column_metadata, pandas.DataFrame):
            if len(column_metadata) != len(self.index):
                raise ValueError('length of column_metadata dataframe does not match the number of columns'.format(key))
            df = column_metadata
        else:
            raise TypeError('column_metadata must be a dictionary or a {} object'.format(pandas.DataFrame.__mro__[0]))
        return df

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
        return self._alignment.chunk_size

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

    def copy(self):
        return self.__class__(
            self.name, self._alignment.copy(),
            chunk_size=self.chunk_size,
            index=self._index.copy(deep=True),
            metadata=deepcopy(self.metadata),
            column_metadata=self._column_metadata.copy(deep=True)
        )

    def set_chunk_size(self, value, copy=False, recasting_func=None, 
                       reset_index=False):
        if reset_index is not True:
            raise ValueError(
                'cannot change chunk size without resetting the current index')
        aln = self
        if copy is True:
            aln = self.copy()
        # Changing chunk size to invalid size will raise an error
        curr_chunk_size = aln.chunk_size
        aln._alignment.chunk_size = value

        # Adjust column metadata
        if recasting_func is None:
            # current size -> 1's -> new size
            df = self._default_expander_func(
                aln._column_metadata, curr_chunk_size)
            print(curr_chunk_size, value, len(df))
            df = self._default_reducer_func(df, value)
            print(curr_chunk_size, value, len(df))
        else:
            df = recasting_func(
                aln._column_metadata, curr_chunk_size, value)
        aln._column_metadata = df
        if copy is True:
            return aln

    def reset_index(self):
        self.cols.reset_index()

    def variants(self):
        return list(self.cols.map(Counter))

    def consensus(self, half=True):
        cons = []
        for cnts in self.cols.map(Counter):
            char, cnt = max(cnts.items(), key=lambda x: x[1])
            if half is True and cnt < self.nrows:
                char = None
            cons.append(char)
        return cons

    def drop(self, value, case_sensitive=False, copy=False, dry_run=False):
        return self.cols.drop(value, case_sensitive=case_sensitive,
                              copy=copy, dry_run=dry_run)

    def drop_except(self, value, case_sensitive=False, copy=False,
                    dry_run=False):
        return self.cols.drop_except(value, case_sensitive=case_sensitive,
                                     copy=copy, dry_run=dry_run)

    def drop_n(self, n_char='N', case_sensitive=False, copy=False,
               dry_run=False):
        return self.cols.drop_n(n_char=n_char, case_sensitive=case_sensitive,
                                copy=copy, dry_run=dry_run)

    def drop_gap(self, gap_char='-', copy=False, dry_run=False):
        return self.cols.drop_gap(gap_char=gap_char, copy=copy, dry_run=dry_run)

    @staticmethod
    def _default_expander_func(df, n):
        return pandas.DataFrame(
                {col: numpy.repeat(df.values, n) for col in df})

    @staticmethod
    def _default_reducer_func(df, n):
        def apply_func(x, col):
            if x[col].dtype != numpy.dtype('O'):
                return max(Counter(x[col]).items(), key=lambda y: y[1])[0]
            return numpy.mean(x[col])
        grouper = numpy.arange(len(df))//n
        return df.groupby(grouper) \
                 .apply(lambda x: pandas.Series(
                    {col: apply_func(x, col) for col in df}))

    # TODO: implement __copy__ and __deepcopy__

    # Special methods
    # ==========================================================================

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.rows.__getitem__(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, str) for val in key)):
            return self.rows.__getitem__(key)
        elif isinstance(key, int) or isinstance(key, slice):
            return self.cols.__getitem__(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, int) for val in key)):
            return self.cols.__getitem__(key)
        raise TypeError('key must be str or int')

    # def __delitem__(self, key):
    #     if isinstance(key, str):
    #         if key in self.ids:
    #             i = self._alignment.row_names_to_indices([key])
    #             return self._alignment.remove_records(i)
    #         raise KeyError('key did not match any identifier')
    #     elif isinstance(key, int):
    #         return self.remove_cols(key)
    #     raise TypeError('key must be str or int')

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
        return hash(self.to_dict(column_metadata=True).items())

    def __eq__(self, other):
        return hash(self) == hash(other)
