from collections import Counter
from copy import copy, deepcopy
import os
import inspect
import warnings

import pandas
import numpy

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import Record
from alignmentrs.util import idseq_to_display
from alignmentrs.aln.mixins import (RecordsSerdeMixin, FastaSerdeMixin,
                                    JsonSerdeMixin, PickleSerdeMixin)
from alignmentrs.history import History
from alignmentrs.history import Record as Record_
from alignmentrs.warning import NoNameWarning, DuplicateNameWarning
from alignmentrs.util import add_to_history
from .mutator import RowData, ColData
from .metadata import MetadataRedirect


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
    def __init__(self, records, name=None, index=None, comments: dict=None,    
                 row_metadata=None, column_metadata=None, store_history=True,
                 **kwargs):
        """Creates a new Alignment object from a sample BaseAlignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        records : list of Record
            Aligned sequences as a list of records.
        index : pandas.Index optional
            Index for alignment columns.
        comments : dict, optional
            Other information related to the alignment. (default is None,
            which creates a blank dictionary)
        row_metadata : pandas.DataFrame or dict, optional
            DataFrame containing annotation for columns.
        column_metadata : pandas.DataFrame or dict, optional
            DataFrame containing annotation for columns.
        store_history : bool, optional
            Whether or not to store actions when the state of the Alignment changes.
        **kwargs
            Other keyword arguments used to initialize states in the
            Alignment object.

        """
        self.name = name
        self.data: BaseAlignment = self.data_constructor(records)
        self.row_metadata = \
            self._row_metadata_constructor(records, row_metadata)
        if index is None:
            index = pandas.Index(range(self.data.ncols))
        self.column_metadata = \
            self._col_metadata_constructor(records, column_metadata, index)
        self.comments = self._comments_constructor(comments)
        self._history = History() if store_history else None
        add_to_history(
            self, self.__class__.__name__, name, records,
            index=index, comments=comments, 
            row_metadata=row_metadata, column_metadata=column_metadata,
            store_history=store_history,
            **kwargs
        )
        self.row = RowData(self)
        self.col = ColData(self)

    # Constructors
    def data_constructor(self, records):
        # Constructs a BaseAlignemnt from the input
        # `records` can be a list of Record or a BaaseAlignment
        # Otherwise, raises a TypeError
        # TODO: Handle possibility of empty alignment
        if isinstance(records, list):
            if sum((isinstance(rec, Record) for rec in records)) \
                == len(records):
                return BaseAlignment(records)
            raise TypeError('records must be a list of Record objects')
        elif isinstance(records, BaseAlignment):
            return records.copy()
        raise TypeError('records must be a list of Record objects or a BaseAlignment')

    def _comments_constructor(self, comments):
        # Constructs metadata dictionary
        # `metadata` can be a ditionary or None
        # TODO: Handle any kind of mapping instead of only a dictionary
        if comments is None:
            return dict()
        elif isinstance(comments, dict):
            return comments
        raise TypeError('comments must be a dictionary of keys and values')

    def _col_metadata_constructor(self, records, column_metadata, index):
        # Constructs column metadata DataFrame from `column_metadata` and
        # `index` inputs`.
        # `column_metadata` must be None, dict, or pandas.DataFrame
        # Otherwise raises TypeError.
        if column_metadata is None:
            df = pandas.DataFrame(None, index=index)
        elif isinstance(column_metadata, dict):
            # Check if values match the length of the index
            for key, val in column_metadata.items():
                if len(val) != len(index):
                    raise ValueError('{} value length does not match the number of columns')
            df = pandas.DataFrame(column_metadata, index=index)
        elif isinstance(column_metadata, pandas.DataFrame):
            if len(column_metadata) != len(index):
                raise ValueError('length of column_metadata dataframe does not match the number of columns')
            df = column_metadata
            df.index = index
        else:
            raise TypeError('column_metadata must be a dictionary or a {} object'.format(pandas.DataFrame.__mro__[0]))
        return df

    def _row_metadata_constructor(self, records, row_metadata):
        # Constructs row metadata DataFrame from `row_metadata`
        index = pandas.Index([rec.id for rec in records])
        df = pandas.DataFrame({
            'description': [rec.description for rec in records]
        }, index=index)
        if isinstance(row_metadata, dict):
            # Check if values match the length of the index
            for key, val in row_metadata.items():
                if len(val) != len(self.nrows):
                    raise ValueError('{} value length does not match the number of rows'.format(key))
            df = df.join(pandas.DataFrame(row_metadata, index=index))
        elif isinstance(row_metadata, pandas.DataFrame):
            if not all(row_metadata.index == index):
                raise ValueError('index of row_metadata DataFrame does not match the ids in the alignment')
            if 'description' in row_metadata.columns:
                df = row_metadata
            else:
                df = df.join(row_metadata)
        elif row_metadata is None:
            pass
        else:
            raise TypeError('row_metadata must be a dictionary or a {} object'.format(pandas.DataFrame.__mro__[0]))
        return df

    # Properties
    @property
    def index(self):
        """pandas.core.indexes.base.Index: Returns the column index
        of the alignment."""
        return self.column_metadata.index

    @index.setter
    def index(self, index):
        """Sets the column index of the alignment."""
        index = pandas.Index(index)
        self.column_metadata.index = index
        # Add history by default
        add_to_history(
            self, '.index', index
        )

    @property
    def nrows(self):
        """int: Returns the number of rows in the alignment."""
        return self.data.nrows

    @property
    def ncols(self):
        """int: Returns the number of columns in the alignment."""
        return self.data.ncols

    @property
    def records(self):
        """list of Record: Returns the list of records."""
        return [Record(vals[0], vals[1]['description'], self.data.get_row(i))
                for i, vals in enumerate(self.row_metadata.iterrows())]

    @property
    def ids(self):
        """list of str: Returns the list of identifiers."""
        return self.row_metadata.index.to_list()

    @property
    def sequences(self):
        """list of str: Returns the list of sequences."""
        return self.data.sequences

    @property
    def history(self):
        """History: Returns history of previous actions performed that may have
        changed the state of the alignment."""
        return self._history

    @property
    def row_and_metadata(self):
        df = self.row_metadata.copy(deep=True)
        df['sequence'] = self.data.sequences
        return df

    @property
    def column_and_metadata(self):
        df = self.column_metadata.copy(deep=True)
        df['sequence'] = [
            ''.join(seq_vec)
            for seq_vec in self.data.get_cols(list(range(self.ncols)))
        ]
        return df


    # Methods
    # ==========================================================================

    def set_record_as_column_metadata(self, i, func, name=None, copy=False,
                                      **kwargs):
        """Transforms a record into column metadata. Removes the record from the
        alignment.
        
        Parameters
        ----------
        i : int
            Index or record identifier.
        func : function
            Function used to transform the record's sequence into
            column metadata values. The functions takes a chunk of sequence as
            input and returns a corresponding value.
        name : str, optional
            Name of the column metadata. (default is None, if not specified,
            the record ID will be used as the name.)
        copy : bool, optional
            Whether to perform the operation on a new copy or 
            do it inplace. (default is False, modifies the alignment inplace)

        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment. Otherwise, the operation is performed inplace and does not return any value.

        """
        aln = self
        if copy is True:
            aln = self.copy()

        # Add to column metadata
        data = [func(v) for v in list(self.data.get_row(i))]        
        if name is None:
            name = self.row[i].id
        self.column_metadata[name] = data
        self.row.remove(i)

        # Add to history
        func_sig = func.__qualname__ + \
            repr(inspect.signature(func)) \
                .lstrip('<Signature ').rstrip('>')
        add_to_history(
            aln, '.set_record_as_column_metadata', i, func_sig,
            name=name,
            copy=copy,
            **kwargs
        )
        # TODO: Remove from row metadata
        if copy is True:
            return aln        

    def copy(self):
        """Creates a deep copy of the alignment.
        
        Returns
        -------
        Alignment
            The new alignment is an new independent copy.
            Changes made in the new alignment will not affect the original, and
            vice versa.

        """
        return deepcopy(self)

    def reset_index(self, copy=False, **kwargs):
        """Resets the alignment index.
        
        Parameters
        ----------
        copy : bool, optional
            Whether to reset the index on a copy of the alignment or
            reset the index inplace. (default is False, the index of the
            alignment is reset inplace)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment with the
            new index. Otherwise, the operation is performed inplace and does
            not return any value.

        """
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        aln.cols.reset_index(copy=False, _record_history=False)
        # Add to history
        add_to_history(aln, '.reset_index', copy=copy, **kwargs)
        if copy is True:
            return aln

    def variants(self):
        """Returns a list of Counter dictionaries containg the identities and
        frequency of variants in each column of the alignment.

        This is a shortcut for using the Counter constructor on the
        `map` method for columns in the alignment.
        
        Returns
        -------
        list of Counter 

        """
        return list(self.col.map(Counter))

    def consensus(self, threshold=0.5):
        """Returns the consensus sequence of the alignment.
        
        Parameters
        ----------
        threshold : bool, optional
            Percent threshold of rows a variant must occur in order to be
            considered the consensus sequence for the column.
            If no variant meets this threshold, `None` is used instead.
            (default is 0.5, which means variants occuring greater than
            or equal to 50% of the rows become the consensus sequence.)
        
        Returns
        -------
        list of str

        """
        cons = []
        for cnts in self.col.map(Counter):
            char, cnt = max(cnts.items(), key=lambda x: x[1])
            if cnt < self.nrows * threshold:
                char = None
            cons.append(char)
        return cons

    def drop(self, value, case_sensitive=False, copy=False, dry_run=False, 
             mode='any', **kwargs):
        """Drops columns given a certain value.
        
        Parameters
        ----------
        value : str
            Character/s to match against.
        case_sensitive : bool, optional
            Whether of not the character search is case sensitive.
            (default is False, search is case insensitive)
        copy : bool, optional
            Whether to drop columns on a copy of the alignment or
            perform the dropping inplace. 
            (default is False, dropping is performed inplace)
        dry_run : bool, optional
            If True, shows the number of columns to be dropped (True)
            and returns the list of column positions to be dropped.
            No changes are made to the alignment at this point.
            If False, columns with the value are removed.
            (default is False, dropping is performed immediately)
        mode : str, optional
            Method used to decide when a column is dropped.
            If the value is 'any', the column will be dropped when
            at least one match is found.
            If the value if 'all', the column will be dropped only if
            the entire column matches the specified value.
            (the default is 'any', any match drops the column)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment after dropping
            matching columns. Otherwise, dropping is performed inplace
            and does not return any value.

        """
        res = self.col.drop(value, case_sensitive=case_sensitive,
                              copy=copy, dry_run=dry_run, mode=mode,
                              _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self, '.drop', value,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def drop_except(self, value, case_sensitive=False, copy=False,
                    dry_run=False, mode='any', **kwargs):
        """Drops columns given a certain value.
        
        Parameters
        ----------
        value : str
            Character/s to match against.
        case_sensitive : bool, optional
            Whether of not the character search is case sensitive.
            (default is False, search is case insensitive)
        copy : bool, optional
            Whether to drop columns on a copy of the alignment or
            perform the dropping inplace. 
            (default is False, dropping is performed inplace)
        dry_run : bool, optional
            If True, shows the number of columns to be kept (True)
            and returns the list of column positions to be kept.
            No changes are made to the alignment at this point.
            If False, columns matching the value are kept.
            (default is False, dropping is performed immediately)
        mode : str, optional
            Method used to decide when a column is dropped.
            If the value is 'any', the column will be kept when
            at least one match is found.
            If the value if 'all', the column will be kept only if
            the entire column matches the specified value.
            (the default is 'any', any match drops the column)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment after dropping
            columns that did not match the value. Otherwise, dropping is
            performed inplace and does not return any value.

        """
        res = self.col.drop_except(value, case_sensitive=case_sensitive,
                                    copy=copy, dry_run=dry_run, mode=mode,
                                    _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self, '.drop_except', value,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def drop_n(self, n_char='N', case_sensitive=False, copy=False,
               dry_run=False, **kwargs):
        """Drops columns that has the 'N' ambiguity character.
        
        Parameters
        ----------
        n_char : str, optional
            Ambiguity character to match against.
            (default is 'N')
        case_sensitive : bool, optional
            Whether of not the character search is case sensitive.
            (default is False, search is case insensitive)
        copy : bool, optional
            Whether to drop columns on a copy of the alignment or
            perform the dropping inplace. 
            (default is False, dropping is performed inplace)
        dry_run : bool, optional
            If True, shows the number of columns to be dropped (True)
            and returns the list of column positions to be dropped.
            No changes are made to the alignment at this point.
            If False, columns with the value are removed.
            (default is False, dropping is performed immediately)
        mode : str, optional
            Method used to decide when a column is dropped.
            If the value is 'any', the column will be dropped when
            at least one match is found.
            If the value if 'all', the column will be dropped only if
            the entire column matches the specified value.
            (the default is 'any', any match drops the column)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment after dropping
            columns with the ambiguity character. Otherwise, dropping is
            performed inplace and does not return any value.

        """
        res = self.col.drop_n(n_char=n_char, case_sensitive=case_sensitive,
                               copy=copy, dry_run=dry_run, 
                               _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self, '.drop_n',
                n_char=n_char,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                **kwargs
            )
        return res

    def drop_gap(self, gap_char='-', copy=False, dry_run=False, **kwargs):
        """Drops columns that has the '-' gap character.
        
        Parameters
        ----------
        n_char : str, optional
            Gap character to match against. (default is '-')
        case_sensitive : bool, optional
            Whether of not the character search is case sensitive.
            (default is False, search is case insensitive)
        copy : bool, optional
            Whether to drop columns on a copy of the alignment or
            perform the dropping inplace. 
            (default is False, dropping is performed inplace)
        dry_run : bool, optional
            If True, shows the number of columns to be dropped (True)
            and returns the list of column positions to be dropped.
            No changes are made to the alignment at this point.
            If False, columns with the value are removed.
            (default is False, dropping is performed immediately)
        mode : str, optional
            Method used to decide when a column is dropped.
            If the value is 'any', the column will be dropped when
            at least one match is found.
            If the value if 'all', the column will be dropped only if
            the entire column matches the specified value.
            (the default is 'any', any match drops the column)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment after dropping
            columns with the gap character. Otherwise, dropping is
            performed inplace and does not return any value.

        """
        res = self.col.drop_gap(gap_char=gap_char, copy=copy, 
                                 dry_run=dry_run, _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self, '.drop_gap',
                gap_char=gap_char,
                copy=copy,
                dry_run=dry_run,
                **kwargs
            )
        return res

    def join(self, others, reset_index=False, copy=False, **kwargs):
        """Marges the current alignment with one or more other alignments.
        This extends the number of columns in the alignment.
        
        Parameters
        ----------
        others : Alignment or list of Alignment
            Other alignments to append to the current alignment.
        reset_index : bool, optional
            Whether to reset the index on a copy of the alignment or
            reset the index inplace. (default is False, the index of the
            alignment is reset inplace)
        copy : bool, optional
            Whether to keep the current alignment intact and create a new copy of the joined alignment or join the alignments inplace.
            (default is False, joining is performed inplace)
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment after joining with the other alignment/s. Otherwise, dropping is
            performed inplace and does not return any value.

        """
        aln = self
        if copy is True:
            aln = self.copy()
        if isinstance(others, Alignment):
            balns = [others.data]
            others = [others]
        elif isinstance(others, list) and \
            sum((isinstance(o, Alignment) for o in others)) == len(others):
            balns = [o.data for o in others]
        else:
            raise ValueError(
                'others must be an Alignment or ''a list of Alignment objects')
        # check if chunks are the same
        # check if number of records are the same
        curr_ncols = aln.ncols
        aln.data.concat(balns)
        # Concat dataframes
        aln.column_metadata = pandas.concat(
            [aln.column_metadata] + 
            [aln.column_metadata for aln in others],
            sort=False, axis=0
        )
        # Add names
        name_list = ([aln.name]*curr_ncols) + \
                    [o.name for o in others for _ in range(o.ncols)]
        if None in name_list:
            warnings.warn('used `None` in _src_name column metadata '
                'because some alignments have no name', NoNameWarning)
        if set(name_list) != len(others) + 1:
            warnings.warn('some alignments have the same name', 
                DuplicateNameWarning)
        aln.column_metadata['_src_name'] = name_list
        # Concat index
        if reset_index is True:
            aln._index = pandas.Index(range(
                sum([len(aln._index)] + [len(o._index) for o in others])
            ))
            aln.column_metadata.reset_index(drop=True, inplace=True)
        else:
            aln._index = pandas.Index(pandas.concat(
                [pandas.Series(aln._index)] + 
                [pandas.Series(o._index) for o in others]
            ))
        # Add to history
        add_to_history(
            aln, '.join', others,
            reset_index=reset_index,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    # Internal methods
    # ==========================================================================

    @staticmethod
    def _default_expander_func(df, n):
        # The default way to expand the dataframe into a dataframe with more
        # rows is to replicate the dataframe by some multiple.
        return pandas.DataFrame(
                {col: numpy.repeat(df.values, n) for col in df})

    @staticmethod
    def _default_reducer_func(df, n):
        # The default way to reduce the dataframe into a dataframe with fewer
        # rows is to get the mean if the data is numerical, otherwise get
        # the most frequently occurring (mode).
        def apply_func(x, col):
            if x[col].dtype != numpy.dtype('O'):
                return max(Counter(x[col]).items(), key=lambda y: y[1])[0]
            return numpy.mean(x[col])
        grouper = numpy.arange(len(df))//n
        return df.groupby(grouper) \
                 .apply(lambda x: pandas.Series(
                    {col: apply_func(x, col) for col in df}))

    def data_state(self):
        return {
            'name': self.name,
            'nrows': self.nrows,
            'ncols': self.ncols,
        }

    def _metadata_state(self):
        return {
            'num_comments': len(self.comments),
            'num_row_metadata': len(self.row_metadata.columns),
            'num_column_metadata': len(self.column_metadata.columns),
        }

    def _state(self):
        return dict(**self.data_state(), **self._metadata_state())

    # Python magic methods
    # ==========================================================================

    def __getitem__(self, key):
        # Allows access to records and columns by indexing
        # If the key is a str or list of str, this is interpreted
        # to mean that records should be returned.
        # If the key is an int or list of ints, this is interpreted
        # that columns should be returned.
        if isinstance(key, str):
            return self.row.__getitem__(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, str) for val in key)):
            return self.row.__getitem__(key)
        elif isinstance(key, int) or isinstance(key, slice):
            return self.col.__getitem__(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, int) for val in key)):
            return self.col.__getitem__(key)
        raise TypeError('key must be str or int')

    # def __delitem__(self, key):
    #     if isinstance(key, str):
    #         if key in self.ids:
    #             i = self.data.row_names_to_indices([key])
    #             return self.data.remove_records(i)
    #         raise KeyError('key did not match any identifier')
    #     elif isinstance(key, int):
    #         return self.remove_cols(key)
    #     raise TypeError('key must be str or int')

    def __iter__(self):
        # Does not implement __iter__ because its meaning is ambiguous.
        raise NotImplementedError(
            'Use .rows.iter() or .rows.iter_sequences() '
            'to iterate over records or sample sequences, respectively.\n'
            'Use .cols.iter() to iterate across columns in the alignment.')

    def __repr__(self):
        # Returns the stringed representation of the alignment.
        parts = []
        parts.append('[Alignment]')
        parts.append('name = {}'.format(self.name))
        parts.append('nrows = {}'.format(self.nrows))
        parts.append('ncols = {}'.format(self.ncols))
        if self:
            aln = idseq_to_display(self.ids, self.sequences)
            parts += ['', aln, '']
        parts.append('[Alignment.Metadata]')
        parts.append('comment_keys = [{}]'.format(
            ', '.join(list(self.comments.keys()))
        ))
        parts.append('row_metadata_keys = [{}]'.format(
            ', '.join(list(self.row_metadata.keys()))
        ))
        parts.append('column_metadata_keys = [{}]'.format(
            ', '.join(list(self.column_metadata.keys()))
        ))
        return '\n'.join(parts)

    def __str__(self):
        # Returns the string representation of the alignment used for printing.
        return str(self.data)

    def __len__(self):
        # len() is not implemented becuase its meaning is ambiguous
        # for an alignment
        raise NotImplementedError(
            'Use .nrows to get the number of samples, or '
            '.ncols to get the number of columns in the alignment.')

    def __bool__(self):
        # Returns True if the instance has one or more records,
        # and that the records have one or more columns
        if self.ncols == 0 or self.nrows == 0:
            return False
        return True

    def __hash__(self):
        # Generates an integer hash using the dictionary representation
        # of the instance
        return hash(self.to_dict(column_metadata=True).items())

    def __eq__(self, other):
        # Implements ability to compare to other objects using a hash
        # May give the same hash value but not necessarily be the same object
        if isinstance(other, self.__class__):
            return hash(self) == hash(other)
        return False

    def __deepcopy__(self, memo):
        # (self, records, name=None, index=None, comments: dict=None,    
        #          row_metadata=None, column_metadata=None, store_history=True,
        #          **kwargs)
        obj = self.__class__.__new__(self.__class__)

        if '_instance' in self.__dict__.keys():
            obj = self._instance.__class__.__new__(self._instance.__class__)
            self = self._instance

        obj.data = self.data.copy()
        obj.name = deepcopy(self.name, memo)
        obj.row_metadata = self.row_metadata.copy(deep=True)
        obj.column_metadata = self.column_metadata.copy(deep=True)
        obj.comments = deepcopy(self.comments, memo)
        obj._history = deepcopy(self._history, memo)
        obj.row = RowData(obj)
        obj.col = ColData(obj)
        
        # obj._history = deepcopy(self.history, memo)
        # TODO: Add history after deepcopy
        # Add to history
        # if obj._history is not None:
        #     obj._history = deepcopy(self._history, memo)
        #     obj._history.add('deepcopy', obj._state())
        return obj
