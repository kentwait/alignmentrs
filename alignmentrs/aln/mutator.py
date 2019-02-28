from copy import deepcopy
import numbers
import inspect

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import Record
from alignmentrs.util import add_to_history


__all__ = ['RowData', 'ColData']


class RowData:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 0

    # def insert(self, position, records, copy=False, **kwargs):
    #     """Inserts one or more records into the alignment naively (without realignment).
        
    #     Parameters
    #     ----------
    #     position : int or str
    #         Position to be inserted at. If str, 
    #     records : Record or list of Record
    #         Record/s to be inserted.
    #     copy : bool, optional
    #         Whether to insert records on a copy of the alignment, keeping
    #         the current alignment unchanged, or insert records inplace. 
    #         (default is False, insertiong is performed inplace)
        
    #     Raises
    #     ------
    #     TypeError
    #         Value of records is not a Record or List of Record
        
    #     Returns
    #     -------
    #     Alignment or None
    #         If copy is True, returns a deep copy of the Alignment with the
    #         inserted records. Otherwise, insertion is performed inplace
    #         and does not return any value.

    #     """
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     if isinstance(records, Record):
    #         aln.data.insert_records(position, records)
    #     elif isinstance(records, list) and \
    #         sum((isinstance(rec, Record) for rec in records)):
    #         aln.data.insert_records(position, records)
    #     else:
    #         raise TypeError('records must be a Record or a list of Record objects')
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.row.insert',
    #         position, records,
    #         copy=copy,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    # def prepend(self, records, copy=False, **kwargs):
    #     """Inserts one or more records before all the existing records in the
    #     alignment naively (without realignment).
        
    #     Parameters
    #     ---------- 
    #     records : Record or list of Record
    #         Record/s to be inserted.
    #     copy : bool, optional
    #         Whether to insert records on a copy of the alignment, keeping
    #         the current alignment unchanged, or insert records inplace. 
    #         (default is False, insertiong is performed inplace)
        
    #     Raises
    #     ------
    #     TypeError
    #         Value of records is not a Record or List of Record
        
    #     Returns
    #     -------
    #     Alignment or None
    #         If copy is True, returns a deep copy of the Alignment with the
    #         inserted records. Otherwise, insertion is performed inplace
    #         and does not return any value.

    #     """
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     aln.insert(0, records, copy=copy, _record_history=False)
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.row.prepend', records,
    #         copy=copy,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    # def append(self, records, copy=False, **kwargs):
    #     """Inserts one or more records after all the existing records in the
    #     alignment naively (without realignment).
        
    #     Parameters
    #     ---------- 
    #     records : Record or list of Record
    #         Record/s to be inserted.
    #     copy : bool, optional
    #         Whether to insert records on a copy of the alignment, keeping
    #         the current alignment unchanged, or insert records inplace. 
    #         (default is False, insertiong is performed inplace)
        
    #     Raises
    #     ------
    #     TypeError
    #         Value of records is not a Record or List of Record
        
    #     Returns
    #     -------
    #     Alignment or None
    #         If copy is True, returns a deep copy of the Alignment with the
    #         inserted records. Otherwise, insertion is performed inplace
    #         and does not return any value.

    #     """
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     if isinstance(records, Record):
    #         aln.data.append_records(records)
    #     elif isinstance(records, list) and \
    #         sum((isinstance(rec, BaseAlignment) for rec in records)):
    #         aln.data.append_record(records)
    #     else:
    #         raise TypeError('records must be a Record or a list of Record objects')
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.row.append', records,
    #         copy=copy,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    def remove(self, positions, copy=False, **kwargs):
        """Removes one or more records from the alignment naively
        (without realignment).
        
        Parameters
        ---------- 
        positions : int, list of int
            Positions to remove.
        copy : bool, optional
            Whether to remove records from a copy of the alignment, keeping
            the current alignment intact, or remove the records inplace. 
            (default is False, insertiong is performed inplace)
        
        Raises
        ------
        TypeError
            Value of records is not a Record or List of Record
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment, removing
            the specified records. Otherwise, removal is performed inplace
            and does not return any value.

        """
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # TODO: Handle str, list of str
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')
        aln.data.remove_rows(positions)
        indices = aln.row_metadata.index[positions]
        aln.row_metadata.drop(indices, axis=0, inplace=True)
        # Add to history
        add_to_history(
            self._instance, '.row.remove', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    def retain(self, positions, copy=False, **kwargs):
        """Retains one or more records from the alignment, and
        removes all the others.
        
        Parameters
        ---------- 
        positions : int, list of int
            Positions to retain.
        copy : bool, optional
            Whether to remove records from a copy of the alignment, keeping
            the current alignment intact, or remove the records inplace. 
            (default is False, insertiong is performed inplace)
        
        Raises
        ------
        TypeError
            Value of records is not a Record or List of Record
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment containing
            only the specified records. Otherwise, removal is performed inplace
            and does not return any value.

        """
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # TODO: Handle str, list of str
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')
        aln.data.retain_rows(positions)
        aln.row_metadata = aln.row_metadata.iloc[positions]
        # Add to history
        add_to_history(
            self._instance, '.row.retain', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    # def drain(self, positions, **kwargs):
    #     """Removes one or more records from the alignment, and returns
    #     these records as a new alignment.
        
    #     Parameters
    #     ---------- 
    #     positions : int, list of int
    #         Positions to drain from the alignment.
        
    #     Raises
    #     ------
    #     TypeError
    #         Value of records is not a Record or List of Record
        
    #     Returns
    #     -------
    #     Alignment
    #         Returns a deep copy of the Alignment containing the drained
    #         records.

    #     """
    #     if isinstance(positions, int):
    #         positions = [positions]
    #     elif isinstance(positions, list) and \
    #         sum((isinstance(pos, int) for pos in positions)) == len(positions):
    #         pass
    #     else:        
    #         raise TypeError('positions must be an int or a list of int')
    #     remove_positions = self._instance.data.invert_rows(positions)
    #     new_baln = self._instance.data.drain_records(remove_positions)
    #     aln = self._instance.__class__(
    #         self._instance.name,
    #         new_baln, 
    #         chunk_size=self._instance.chunk_size,
    #         index=self._instance._index.copy(deep=True), 
    #         metadata=deepcopy(self._instance.metadata), 
    #         column_metadata=self._instance.column_metadata.copy(deep=True))
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.row.drain[from]', positions,
    #         **kwargs
    #     )
    #     add_to_history(
    #         aln, '.row.drain[to]', positions,
    #         **kwargs
    #     )
    #     return aln

    def reorder(self, positions, copy=False, **kwargs):
        """Reorder records in the alignment based on a given order.
        
        Parameters
        ---------- 
        positions : list of int
            Position list.
        copy : bool, optional
            Whether to remove records from a copy of the alignment, keeping
            the current alignment intact, or remove the records inplace. 
            (default is False, insertiong is performed inplace)
        
        Raises
        ------
        TypeError
            Value of records is not a Record or List of Record
        
        Returns
        -------
        Alignment or None
            If copy is True, returns a deep copy of the Alignment containing
            only the specified records. Otherwise, removal is performed inplace
            and does not return any value.

        """
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:        
            raise TypeError('positions must be an int or a list of int')
        aln.data.reorder_rows(positions)
        aln.row_metadata = aln.row_metadata.iloc[positions]
        # Add to history
        add_to_history(
            self._instance, '.row.reorder', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    def filter(self, function, copy=False, dry_run=False, inverse=False,
               **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # Function accepts a Record, outputs true or false
        if not(function is not None and callable(function)):
            raise TypeError('missing filter function')
        positions = [i for i in range(aln.nrows) 
                     if function(aln.data.get_row(i))]
        remove_positions = aln.data.invert_rows(positions)
        if dry_run:
            parts = []
            parts.append('[Filter]')
            parts.append('True = {}/{}'.format(
                len(positions), aln.nrows))
            parts.append('False = {}/{}'.format(
                len(remove_positions), aln.nrows))
            print('\n'.join(parts))
            return {
                'function': function,
                True: positions,
                False: remove_positions
            }
        if inverse:
            aln.row.remove(positions, _record_history=False)
        else:
            aln.row.remove(remove_positions, _record_history=False)
        # Add to history
        func_sig = function.__qualname__ + \
            repr(inspect.signature(function)) \
                .lstrip('<Signature ').rstrip('>')
        add_to_history(
            aln, '.row.filter', func_sig,
            copy=copy,
            dry_run=dry_run,
            inverse=inverse,
            **kwargs
        )
        if copy is True:
            return aln

    def map(self, function):
        for record in self.iter():
            yield function(record)

    def iter(self):
        for i in range(self._instance.nrows):
            yield self._instance.data.get_row(i)

    def iter_sequences(self):
        for i in range(self._instance.nrows):
            yield self._instance.data.get_row(i)

    def __iter__(self):
        return self.iter()

    def __getitem__(self, key):
        # if isinstance(key, str):
        #     if key in self._instance.ids:
        #         i = self._instance.data.row_names_to_indices([key])
        #         return self._instance.data.get_rows([i[0]])[0]
        #     raise KeyError('key did not match any identifier')
        # elif isinstance(key, list) and \
        #     sum((isinstance(val, str) for val in key)):
        #     return self._instance.data.get_records_by_name(key)
        if isinstance(key, int):
            return self._instance.data.get_row(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, int) for val in key)):
            keys = key
        elif isinstance(key, slice):
            abs_slice = key.indices(self._instance.nrows)
            keys = list(range(*abs_slice))
        else:
            raise TypeError('key must be int, list of int, or a slice')
        return self._instance.data.get_rows(keys)

    def __len__(self):
        return self._instance.nrows

    def __repr__(self):
        return self._instance.__repr__()

    def __str__(self):
        return self._instance.__str__()


class ColData:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 1

    # def insert(self, position, values, copy=False,
    #            column_values=None, reset_index=False, **kwargs):
    #     if reset_index is not True:
    #         raise ValueError(
    #             'cannot perform an insert without resetting the current index')
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     # TODO: Check data type of position
    #     if isinstance(values, list) and \
    #         sum((isinstance(val, str) for val in values)):
    #         values = [values]
    #     elif isinstance(values, list) and \
    #         sum((isinstance(val, BaseAlignment)
    #             for lst in values for val in lst)):
    #         pass
    #     else:
    #         raise TypeError(
    #             'values must be a list of str or a list of list of str')
    #     aln.data.insert_cols(position, values)
    #     aln._index, aln.column_metadata = \
    #         self._insert_metadata(aln, position, column_values)
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.col.insert', position, values,
    #         copy=copy,
    #         column_values=column_values,
    #         reset_index=reset_index,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    # def prepend(self, values, copy=False, column_values=None, reset_index=False,
    #             **kwargs):
    #     if reset_index is not True:
    #         raise ValueError(
    #             'cannot perform an prepend without resetting the current index')
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     self.insert(0, values, copy=copy, column_values=column_values,       
    #                 reset_index=reset_index, )
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.col.prepend', values,
    #         copy=copy,
    #         column_values=column_values,
    #         reset_index=reset_index,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    # def append(self, values, copy=False, column_values=None, reset_index=False, 
    #            **kwargs):
    #     if reset_index is not True:
    #         raise ValueError(
    #             'cannot perform an append without resetting the current index')
    #     aln = self._instance
    #     if copy is True:
    #         aln = self._instance.copy()
    #     position = len(aln._index)
    #     if isinstance(values, list) and \
    #         sum((isinstance(val, str) for val in values)):
    #         aln.data.apend_col(values)
    #     elif isinstance(values, list) and \
    #         sum((isinstance(val, BaseAlignment)
    #             for lst in values for val in lst)):
    #         aln.data.apend_cols(values)
    #     else:
    #         raise TypeError(
    #             'values must be a list of str or a list of list of str')
    #     aln.data.insert_cols(position, values)
    #     aln._index, aln.column_metadata = \
    #         self._insert_metadata(aln, position, column_values)
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.col.append', values,
    #         copy=copy,
    #         column_values=column_values,
    #         reset_index=reset_index,
    #         **kwargs
    #     )
    #     if copy is True:
    #         return aln

    def remove(self, positions, copy=False, **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            len(positions) == 0:
            pass
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')
        retain_positions = aln.data.invert_cols(positions)        
        aln.data.remove_cols(positions)
        indices = aln.row_metadata.index[positions]
        aln.column_metadata.drop(indices, axis=0, inplace=True)
        # Add to history
        add_to_history(
            self._instance, '.col.remove', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    def retain(self, positions, copy=False, **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:        
            raise TypeError('positions must be an int or a list of int')
        remove_positions = aln.data.invert_cols(positions)
        aln.data.remove_cols(remove_positions)
        aln.column_metadata = aln.column_metadata.iloc[positions]
        # Add to history
        add_to_history(
            self._instance, '.col.retain', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    # def drain(self, positions, **kwargs):
    #     if isinstance(positions, int):
    #         positions = [positions]
    #     elif isinstance(positions, list) and \
    #         sum((isinstance(pos, int) for pos in positions)):
    #         pass
    #     else:        
    #         raise TypeError('positions must be an int or a list of int')
    #     remove_positions = self._instance.data.invert_cols(positions)
    #     new_baln = self._instance.data.drain_cols(remove_positions)
    #     new_col_metadata = self._instance \
    #         .column_metadata.iloc[remove_positions].copy(deep=True)
    #     self._instance.column_metadata = self._instance \
    #         .column_metadata.iloc[positions].copy(deep=True)
    #     aln = self._instance.__class__(
    #         self._instance.name,
    #         new_baln, 
    #         chunk_size=self._instance.chunk_size,
    #         index=pandas.Index(new_col_metadata.index), 
    #         metadata=deepcopy(self._instance.metadata), 
    #         column_metadata=new_col_metadata)
    #     # Add to history
    #     add_to_history(
    #         self._instance, '.col.drain[from]', positions,
    #         **kwargs
    #     )
    #     add_to_history(
    #         aln, '.col.drain[to]', positions,
    #         **kwargs
    #     )
    #     return aln

    def reorder(self, positions, copy=False, **kwargs):
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
        aln.data.reorder_cols(positions)
        aln.column_metadata = aln.column_metadata.iloc[positions]
        # Add to history
        add_to_history(
            aln, '.col.reorder', positions,
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    def filter(self, function, copy=False, dry_run=False, inverse=False,
               **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        # Function accepts a list of str, outputs true or false
        if not(function is not None and callable(function)):
            raise TypeError('missing filter function')
        positions = [i for i in range(aln.ncols) 
                     if function(aln.data.get_col(i))]
        remove_positions = aln.data.invert_cols(positions)
        if dry_run:
            parts = []
            parts.append('[Filter]')
            parts.append('True = {}/{}'.format(
                len(positions), aln.ncols))
            parts.append('False = {}/{}'.format(
                len(remove_positions), aln.ncols))
            print('\n'.join(parts))
            return {
                'function': function,
                True: positions,
                False: remove_positions
            }
        if inverse:
            aln.col.remove(positions, _record_history=False)
        else:
            aln.col.remove(remove_positions, _record_history=False)
        # Add to history
        func_sig = function.__qualname__ + \
            repr(inspect.signature(function)) \
                .lstrip('<Signature ').rstrip('>')
        add_to_history(
            aln, '.col.filter', func_sig,
            copy=copy,
            dry_run=dry_run,
            inverse=inverse,
            **kwargs
        )
        if copy is True:
            return aln

    def drop(self, value, case_sensitive=False, copy=False, dry_run=False,
             mode='any', **kwargs):
        if len(value) < self._instance.chunk_size:
            if case_sensitive and mode == 'any':
                func = lambda x: \
                    sum([value in chars for chars in x]) > 0
            elif case_sensitive and mode == 'all':
                func = lambda x: \
                    sum([value in chars for chars in x]) == len(x)
            elif not case_sensitive and mode == 'any':
                func = lambda x: \
                    sum([value.upper() in chars.upper() for chars in x]) > 0
            elif not case_sensitive and mode == 'all':
                func = lambda x: \
                    sum([value.upper() in chars.upper() for chars in x]) \
                        == len(x)
            else:
                raise ValueError('invalid mode')
        elif len(value) == self._instance.chunk_size:
            if case_sensitive and mode == 'any':
                func = lambda x: value in x
            elif case_sensitive and mode == 'all':
                func = lambda x: [value]*len(x) == x
            elif not case_sensitive and mode == 'any':
                func = lambda x: \
                    value.upper() in [chars.upper() for chars in x]
            elif not case_sensitive and mode == 'all':
                func = lambda x: \
                    [value.upper()]*len(x) == [chars.upper() for chars in x]
            else:
                raise ValueError('invalid mode')
        else:
            raise ValueError('value is larger than chunk size')
        res = self.filter(func, copy=copy, dry_run=dry_run, inverse=True,
                          _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self._instance, '.col.drop', value,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def drop_except(self, value, case_sensitive=False, copy=False, 
                    dry_run=False, mode='any', **kwargs):
        if len(value) < self._instance.chunk_size:
            if case_sensitive and mode == 'any':
                func = lambda x: \
                    sum([value in chars for chars in x]) > 0
            elif case_sensitive and mode == 'all':
                func = lambda x: \
                    sum([value in chars for chars in x]) == len(x)
            elif not case_sensitive and mode == 'any':
                func = lambda x: \
                    sum([value.upper() in chars.upper() for chars in x]) > 0
            elif not case_sensitive and mode == 'all':
                func = lambda x: \
                    sum([value.upper() in chars.upper() for chars in x]) \
                        == len(x)
            else:
                raise ValueError('invalid mode')
        elif len(value) == self._instance.chunk_size:
            if case_sensitive and mode == 'any':
                func = lambda x: value in x
            elif case_sensitive and mode == 'all':
                func = lambda x: [value]*len(x) == x
            elif not case_sensitive and mode == 'any':
                func = lambda x: \
                    value.upper() in [chars.upper() for chars in x]
            elif not case_sensitive and mode == 'all':
                func = lambda x: \
                    [value.upper()]*len(x) == [chars.upper() for chars in x]
            else:
                raise ValueError('invalid mode')
        else:
            raise ValueError('value is larger than chunk size')
        res = self.filter(func, copy=copy, dry_run=dry_run, inverse=False,
                          _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self._instance, '.col.drop_except', value,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def drop_n(self, n_char='N', case_sensitive=False, copy=False, 
               dry_run=False, mode='any', **kwargs):
        res = self.drop(n_char, case_sensitive=case_sensitive,
                         copy=copy, dry_run=dry_run, mode=mode,
                         _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self._instance, '.col.drop_n',
                n_char=n_char,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def drop_gap(self, gap_char='-', copy=False, dry_run=False, mode='any',
                 **kwargs):
        case_sensitive = True
        res = self.drop(gap_char, case_sensitive=case_sensitive, copy=copy,
                        dry_run=dry_run, mode=mode,
                        _record_history=False)
        # Add to history
        if not dry_run:
            add_to_history(
                res if copy else self._instance, '.col.drop_gap',
                gap_char=gap_char,
                case_sensitive=case_sensitive,
                copy=copy,
                dry_run=dry_run,
                mode=mode,
                **kwargs
            )
        return res

    def map(self, function, skip_n=None, chunk_size=None):
        for col in self.iter(skip_n=skip_n, chunk_size=chunk_size):
            yield function(col)

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
                yield self._instance.data.get_col(i)
            else:
                yield self._instance.data.get_chunk(i, chunk_size)

    def reset_index(self, copy=False, **kwargs):
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        aln.column_metadata.reset_index(drop=True, inplace=True)
        # Add to history
        add_to_history(
            aln, '.col.reset_index',
            copy=copy,
            **kwargs
        )
        if copy is True:
            return aln

    @staticmethod
    def _insert_metadata(aln, position, column_values):
        if isinstance(column_values, list) and \
            sum((isinstance(val, list) for val in column_values)):
            df = pandas.DataFrame(
                {k:v for k,v in zip(aln.column_metadata, column_values)})
        elif isinstance(column_values, list) and \
            sum((isinstance(val, dict) for val in column_values)):
            df = pandas.DataFrame(column_values)
        elif isinstance(column_values, list) and \
            sum((isinstance(val, numbers.Number) or isinstance(val, str)
                for val in column_values)):
            df = pandas.DataFrame(
                {k:v for k,v in zip(aln.column_metadata, column_values)})
        elif isinstance(column_values, dict) and \
            sum((isinstance(val, list) for val in column_values)):
            df = pandas.DataFrame(column_values)
        elif isinstance(column_values, dict) and \
            sum((isinstance(val, numbers.Number) or isinstance(val, str)
                for val in column_values)):
            df = pandas.DataFrame(column_values)
        new_column_metadata = \
            pandas.concat([aln.column_metadata.iloc[:position],
                            df,
                            aln.column_metadata.iloc[position+len(df):]]) \
                .reset_index(drop=True)
        new_index = pandas.Index(len(aln._index) + len(df))
        return new_index, new_column_metadata

    def __iter__(self):
        return self.iter(skip_n=1, chunk_size=1)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._instance.data.get_col(key)
        elif isinstance(key, list) and \
            sum((isinstance(val, int) for val in key)):
            keys = key
        elif isinstance(key, slice):
            abs_slice = key.indices(self._instance.ncols)
            keys = list(range(*abs_slice))
        else:
            raise TypeError('key must be int, list of int, or a slice')
        return self._instance.data.get_cols(keys)

    def __len__(self):
        return self._instance.ncols

    def __repr__(self):
        return self._instance.__repr__()

    def __str__(self):
        return self._instance.__str__()
