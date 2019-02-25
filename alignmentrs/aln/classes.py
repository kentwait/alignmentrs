from collections import Counter
from copy import deepcopy
import warnings
import os
from datetime import datetime
import json
import inspect

import pandas
import numpy

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import BaseRecord
from alignmentrs.util import idseq_to_display
from alignmentrs.aln.mixins import (RecordsSerdeMixin, FastaSerdeMixin,
                                    JsonSerdeMixin, PickleSerdeMixin)
from .mutator import RowMutator, ColMutator


__all__ = ['Alignment', 'CatAlignment']


class NoNameWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass

class DuplicateNameWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass

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
                 store_history=True, **kwargs):
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
        self._history = History() if store_history else None
        # self._store_state_history = store_state_history


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
        elif isinstance(index, list):
            return pandas.Index(index)
        raise TypeError(
            'index must be a list or {} object'.format(pandas.Index.__mro__[0]))

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
            column_metadata.index = self.index
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

    @index.setter
    def set_index(self, index):
        index = pandas.Index(index)
        self._column_metadata.index = index
        self._index = index
        # Add to history
        if self._history is not None:
            self._history.add('.cols.set_index',
                args=[index]
                )

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
    def records(self):
        return self._alignment.records

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

    @property
    def history(self):
        return self._history


    # Methods
    # ==========================================================================

    def copy(self, **kwargs):
        # TODO: Copy over previous history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.copy')
        return self.__class__(
            self.name, self._alignment.copy(),
            chunk_size=self.chunk_size,
            index=self._index.copy(deep=True),
            metadata=deepcopy(self.metadata),
            column_metadata=self._column_metadata.copy(deep=True)
        )

    def set_chunk_size(self, value, copy=False, recasting_func=None, 
                       reset_index=False, **kwargs):
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
        # Add to history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (aln._history is not None):
            aln._history.add('.set_chunk_size',
                args=[value],
                kwargs={
                    'copy': copy,
                    'recasting_func': 
                        repr(inspect.signature(recasting_func))
                            .lstrip('<Signature ').rstrip('>'),
                    'reset_index': reset_index
                })
        if copy is True:
            return aln

    def add_column_metadata(self, name, data, **kwargs):
        if name in self._column_metadata:
            raise ValueError('name already exists')
        self._column_metadata[name] = data
        # Add to history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.set_chunk_size',
                args=[name, data])

    def remove_column_metadata(self, name, **kwargs):
        self._column_metadata.drop(name, axis=1, inplace=True, _record_history=False)
        # Add to history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.remove_column_metadata',
                args=[name])

    def reset_index(self, copy=False, **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # Add to history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            aln._history.add('.reset_index')
        aln.cols.reset_index(copy=False, _record_history=False)
        if copy is True:
            return aln

    def variants(self):
        return list(self.cols.map(Counter))

    def consensus(self, half=True):
        cons = []
        for cnts in self.cols.map(Counter):
            char, cnt = max(cnts.items(), key=lambda x: x[1])
            if half is True and cnt < self.nrows/2:
                char = None
            cons.append(char)
        return cons

    def drop(self, value, case_sensitive=False, copy=False, dry_run=False, 
             mode='any', **kwargs):
        # Add to history
        record = True
        if dry_run:
            record = False
        elif '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.drop',
                args=[value],
                kwargs={
                    'case_sensitive': case_sensitive,
                    'copy': copy,
                    'dry_run': dry_run,
                })
        return self.cols.drop(value, case_sensitive=case_sensitive,
                              copy=copy, dry_run=dry_run, mode=mode,
                              _record_history=False)

    def drop_except(self, value, case_sensitive=False, copy=False,
                    dry_run=False, mode='any', **kwargs):
        # Add to history
        record = True
        if dry_run:
            record = False
        elif '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.drop_except',
                args=[value],
                kwargs={
                    'case_sensitive': case_sensitive,
                    'copy': copy,
                    'dry_run': dry_run,
                    'mode': mode,
                })
        return self.cols.drop_except(value, case_sensitive=case_sensitive,
                                     copy=copy, dry_run=dry_run, mode=mode,
                                     _record_history=False)

    def drop_n(self, n_char='N', case_sensitive=False, copy=False,
               dry_run=False, **kwargs):
        # Add to history
        record = True
        if dry_run:
            record = False
        elif '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.drop_n',
                kwargs={
                    'n_char': n_char,
                    'case_sensitive': case_sensitive,
                    'copy': copy,
                    'dry_run': dry_run,
                })
        return self.cols.drop_n(n_char=n_char, case_sensitive=case_sensitive,
                                copy=copy, dry_run=dry_run, _record_history=False)

    def drop_gap(self, gap_char='-', copy=False, dry_run=False, **kwargs):
        # Add to history
        record = True
        if dry_run:
            record = False
        elif '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (self._history is not None):
            self._history.add('.drop_gap',
                kwargs={
                    'gap_char': gap_char,
                    'copy': copy,
                    'dry_run': dry_run,
                })
        return self.cols.drop_gap(gap_char=gap_char, copy=copy, 
                                  dry_run=dry_run, _record_history=False)

    def join(self, others, reset_index=False, copy=False, **kwargs):
        aln = self
        if copy is True:
            aln = self.copy()
        if isinstance(others, Alignment):
            balns = [others._alignment]
            others = [others]
        elif isinstance(others, list) and \
            sum((isinstance(o, Alignment) for o in others)) == len(others):
            balns = [o._alignment for o in others]
        else:
            raise ValueError(
                'others must be an Alignment or ''a list of Alignment objects')
        # check if chunks are the same
        # check if number of records are the same
        curr_ncols = aln.ncols
        aln._alignment.concat(balns)
        # Concat dataframes
        aln._column_metadata = pandas.concat(
            [aln._column_metadata] + 
            [aln._column_metadata for aln in others],
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
        aln._column_metadata['_src_name'] = name_list
        # Concat index
        if reset_index is True:
            aln._index = pandas.Index(range(
                sum([len(aln._index)] + [len(o._index) for o in others])
            ))
            aln._column_metadata.reset_index(drop=True, inplace=True)
        else:
            aln._index = pandas.Index(pandas.concat(
                [pandas.Series(aln._index)] + 
                [pandas.Series(o._index) for o in others]
            ))
        # Add to history
        record = True
        if '_record_history' in kwargs.keys():
            record = kwargs['_record_history']
            del kwargs['_record_history']
        if record and (aln._history is not None):
            aln._history.add('.join',
                args=[others],
                kwargs={
                    'reset_index': reset_index,
                    'copy': copy,
                })
        if copy is True:
            return aln

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

class History:
    def __init__(self, *args, **kwargs):
        self.records = []
        # self.store_state_history = store_state_history

    def add(self, operation, args=None, kwargs=None,
            # before_state=None, after_state=None
            ):
        # if self.store_state:
        #     if before_state:
        #         raise ValueError('store_state_history is True but "before" state is not specfied')
        #     if after_state:
        #         raise ValueError('store_state_history is True but "after" state is not specfied')
        record = Record(operation, args=args, kwargs=kwargs,
                        # before_state=before_state, after_state=after_state
                        )
        self.records.append(record)
    
    def to_json(self, path=None):
        return '[{}]'.format(
            ', '.join([i.to_json() for i in self.records])
        )

    def to_markdown(self, path=None):
        return '\n'.join([i.to_markdown() for i in self.records])

    def to_string(self, str_formatter=None, strftime='%m/%d/%Y %H:%M:%S'):
        return '\n'.join([
            i.to_string(str_formatter=str_formatter, strftime=strftime)
            for i in self.records
        ])

    def to_toml(self, path=None):
        return '\n'.join([i.to_toml() for i in self.records])

    def to_csv(self, path=None, sep='\t'):
        return '\n'.join([i.to_csv(sep=sep) for i in self.records])

    def __str__(self):
        return '\n'.join([str(item) for item in self.records])

    def __repr__(self):
        return repr(self.records)

    def __getitem__(self, k):
        return self.record[k]

    def __setitem__(self, k, v):
        self.records[k] = v


class Record:
    def __init__(self, operation, args=None, kwargs=None,
                #  before_state=None, after_state=None, 
                 str_formatter=None, module=None):
        self.operation = operation
        self.optype = 'operation'
        if module:
            self.operation = '.'.join([module, self.operation])
        self.args = [repr(arg).replace('\n', ' ') for arg in args] \
            if args else []
        self.kwargs = {key: repr(kwarg).replace('\n', ' ') 
                       for key, kwarg in kwargs.items()} \
            if kwargs else dict()
        self.datetime = datetime.now()
        self._threshold_args = 3
        self._threshold_kwargs = 2
        self.str_formatter = str_formatter

    def to_json(self, path=None):
        return json.dumps(
            {self.optype: 
                {
                    'datetime': self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                    'operation': self.operation,
                    'args': self.args,
                    'kwargs': self.kwargs,
                }
            }
        )

    def to_markdown(self, path=None):
        params = ''
        if len(self.args) > 0:
            params += ', '.join(self.args)
        if len(self.kwargs) > 0:
            if len(self.args) > 0:
                params += ', '
            params += ', '.join(
                ['{}={}'.format(k, v) for k, v in self.kwargs.items()]
            )
        return ('## {optype} {op}\n'
                '  * date and time - {dt}\n'
                '  * statement - `{op}({params})`\n'.format(
                   dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                   op=self.operation,
                   optype=self.optype,
                   params=params,
                ))

    def to_string(self, str_formatter=None, strftime='%m/%d/%Y %H:%M:%S'):
        if str_formatter:
            return str_formatter(
                operation=self.operation,
                optype=self.optype,
                args=self.args,
                kwargs=self.kwargs,
                datetime=self.datetime.strftime(strftime),
            )
        if self.str_formatter:
            return self.str_formatter(
                operation=self.operation,
                optype=self.optype,
                args=self.args,
                kwargs=self.kwargs,
                datetime=self.datetime,
            )
        args = [v for i, v in enumerate(self.args)
                if i < self._threshold_args]
        if len(self.args) > self._threshold_args:
            args += ['...']
        kwargs = ['{k}={v}'.format(k=kv[0], v=kv[1]) 
                for i, kv in enumerate(self.kwargs.items())
                if i < self._threshold_args]
        if len(self.kwargs) > self._threshold_kwargs:
            kwargs += ['...']
        params = ''
        if len(args) > 0:
            params += ', '.join(args)
        if len(kwargs) > 0:
            if len(self.args) > 0:
                params += ', '
            params += ', '.join(kwargs)
        return '{optype} {dt} {op}({params})'.format(
            dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
            optype=self.optype,
            op=self.operation,
            params=params,
        )

    def to_toml(self, path=None):
        return ('[{optype}]\n'
                'datetime = {dt}\n'
                'operation = \'{op}\'\n'
                'args = \'{args}\'\n'
                'kwargs = \'{kwargs}\'\n'.format(
                   dt=self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
                   op=self.operation,
                   optype=self.optype,
                   args=self.args,
                   kwargs=str(self.kwargs),
               ))

    def to_csv(self, path=None, sep='\t'):
        return sep.join([
            self.optype,
            self.datetime.strftime('%m/%d/%Y %H:%M:%S'),
            self.operation,
            str(self.args),
            str(self.kwargs),
        ])

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return '<{clsname} datetime={dt} statement={op}({params})>'.format(
            clsname=self.__class__.__name__,
            dt=self.datetime.strftime('%m/%d/%YT%H:%M:%S'),
            op=self.operation,
            params='...' if len(self.args) > 0 or len(self.kwargs) > 0 else ''
        )