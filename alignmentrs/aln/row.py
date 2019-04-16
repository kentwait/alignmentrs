from copy import deepcopy
import numbers
import inspect
import itertools

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import Record
from alignmentrs.utils import add_to_history


__all__ = ['RowData']


class RowMethods:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 0

    def get(self, positions, **kwargs):
        """Returns one or more rows from the alignment as a new alignment.
        
        Parameters
        ----------
        positions : int or iterable
            Integer index positions of rows to return.
        
        Returns
        -------
        Alignment
            Returns the subset of the alignment containing only the specified
            rows. This returns a copy of the original alignment.

        """
        # Check input
        if isinstance(positions, int):
            positions = [positions]
        elif isinstance(positions, list) and \
            sum((isinstance(pos, int) for pos in positions)) == len(positions):
            pass
        else:
            raise TypeError('positions must be an int or a list of int')

        aln = self._instance
        return aln.row.retain(positions, copy=True)

    def remove(self, positions, copy=False, **kwargs):
        """Removes the specified row/s from the alignment naively
        (without realignment).

        Parameters
        ---------- 
        positions : int, list of int
            Integer index positions of rows to remove.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing it done inplace)

        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing the specified rows.

        """
        # Check input
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

        aln = self._instance
        if copy is True:
            aln = self._instance.copy()

        # Remove rows from SeqMatrix
        aln.data.remove_rows(positions)
        # Remove row metadata
        indices = aln.row_metadata.index[positions]
        aln.row_metadata.drop(indices, axis=0, inplace=True)

        # # Add to history
        # add_to_history(
        #     self._instance, '.row.remove', positions,
        #     copy=copy,
        #     **kwargs
        # )

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

