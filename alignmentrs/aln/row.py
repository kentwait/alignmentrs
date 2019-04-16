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
        """Returns the list of row positions where the given function is True.

        Parameters
        ----------
        function : callable
            Function used to evaluate each row. The function should expect a list of str as input and a bool and returns a bool as output.
        copy : bool, optional
            Whether to return a new copy of the filtered alignment, keeping the
            original intact, or filter the alignment inplace.
            (default is False, filtering is performed inplace)
        dry_run : bool, optional
            If True, evaluates the function and returns the list of True and False column position only. Nothing is edited. Otherwise,
            column positions that evaluated False are removed from the
            alignment. (the default is False, the alignment is edited)
        inverse : bool, optional
            If True, columns that evaluate True are removed. (default is False, column positions that evaluate False are removed)

        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing columns that evaluated False.
        dict
            When `dry_run` is True, return a dictionary where True and False
            are keys, and the list of respective column positions are the
            values.

        """
        # Check input

        # Checks if function is callable
        # Function accepts a list of str, outputs true or false
        if not(function is not None and callable(function)):
            raise TypeError('missing filter function')

        # Check optional kwargs
        custom_title = 'Filter'
        if 'custom_title' in kwargs.keys():
            custom_title = kwargs['custom_title']
        custom_class_true = True
        if 'custom_class_true' in kwargs.keys():
            custom_class_true = kwargs['custom_class_true'] + ' (True)'
        custom_class_false = False
        if 'custom_class_false' in kwargs.keys():
            custom_class_false = kwargs['custom_class_false']  + ' (False)'

        aln = self._instance
        if copy is True:
            aln = self._instance.copy()

        # Get the list of positions based on the result of the filtering
        # function. Only position of rows that are True are recorded. 
        positions = [i for i, row in enumerate(aln.data) 
                     if function(row)]
        other_positions = aln.data.invert_rows(positions)

        # `dry_run` shows True or False rows as evaluated by the filtering
        # function.
        # If dry_run=True, this prints out the number of rows in True or
        # False categories and returns the lists of integer positons in
        # True or False classifications.
        if dry_run:
            parts = []
            parts.append('[{}]'.format(custom_title))
            parts.append('{} = {}/{}'.format(
                custom_class_true, len(positions), aln.nrows))
            parts.append('{} = {}/{}'.format(
                custom_class_false, len(other_positions), aln.nrows))
            print('\n'.join(parts))
            return {
                True: positions,
                False: other_positions
            }

        # The default behavior of filter is to keep rows evaluating True for the
        # given filter function, and will remove rows that are False.
        # However, if `inverse` is True, the filter method will do the opposite.
        # It will keep rows that are False and will remove rows that are True.
        if inverse:
            aln.row.remove(positions, _record_history=False)
        else:
            aln.row.retain(positions, _record_history=False)

        # # Add to history
        # func_sig = function.__qualname__ + \
        #     repr(inspect.signature(function)) \
        #         .lstrip('<Signature ').rstrip('>')
        # add_to_history(
        #     aln, '.row.filter', func_sig,
        #     copy=copy,
        #     dry_run=dry_run,
        #     inverse=inverse,
        #     **kwargs
        # )

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

