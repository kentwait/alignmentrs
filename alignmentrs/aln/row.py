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
            Position index/indices of rows to return.
        
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
            Position index/indices of rows to remove.
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
        """Retains the specified row/s in the alignment. Removes all other
        rows naively (without realignment).
        
        Parameters
        ---------- 
        positions : int, list of int
            Position index/indices of rows to retained.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing rows.

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
        if copy is True:
            aln = self._instance.copy()

        # Remove/retain rows
        aln.data.retain_rows(positions)
        aln.row_metadata = aln.row_metadata.iloc[positions]

        # # Add to history
        # add_to_history(
        #     self._instance, '.row.retain', positions,
        #     copy=copy,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def reorder(self, position_list, copy=False, **kwargs):
        """Reorder samples according to the specified list of positions.
        
        Parameters
        ---------- 
        positions : list of int
            Ordered list of integer positons indicating the new order of rows
            in the alignment.
        copy : bool, optional
            Whether to return a new copy of the reordered alignment, keeping the
            original intact, or reorder the alignment inplace. (default is
            False, reordering is done inplace)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after reordering
            rows.

        """
        # Check input
        if isinstance(position_list, list) and \
            sum((isinstance(pos, int) for pos in position_list)):
            if len(list) != self._instance.ncols:
                raise TypeError('length of position list must be equal to the '
                    'number of rows in the alignment: {} != {}'.format(
                        len(list), self._instance.nrows
                    ))
        else:
            raise TypeError('position list must be a list of int')

        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        aln.data.reorder_rows(position_list)
        aln.row_metadata = aln.row_metadata.iloc[position_list]

        # # Add to history
        # add_to_history(
        #     self._instance, '.row.reorder', positions,
        #     copy=copy,
        #     **kwargs
        # )

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

