from copy import deepcopy
import numbers
import inspect
import itertools

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import Record
from alignmentrs.utils import add_to_history


__all__ = ['ColData']


class ColMethods:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 1

    def get(self, positions, **kwargs):
        """Returns one or more columns from the alignment as a new alignment.
        
        Parameters
        ----------
        positions : int or iterable
            Position index/indices of columns to return.
        
        Returns
        -------
        Alignment
            Returns the subset of the alignment containing only
            the specified columns. This returns a copy of the
            original alignment.

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
        return aln.col.retain(positions, copy=True)

    def remove(self, positions, copy=False, **kwargs):
        """Removes the specified column/s from the alignment.
        
        Parameters
        ----------
        positions : int or iterable
            Position index/indices of columns to remove.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing the specified
            columns.

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

        # Remove columns from SeqMatrix
        aln.data.remove_cols(positions)
        # Remove column metadata
        indices = aln.column_metadata.index[positions]
        aln.column_metadata.drop(indices, axis=0, inplace=True)

        # # Add to history
        # add_to_history(
        #     self._instance, '.col.remove', positions,
        #     copy=copy,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def retain(self, positions, copy=False, **kwargs):
        """Retains the specified column/s in the alignment. Removes all the
        other columns.
        
        Parameters
        ----------
        positions : int or iterable
            Position index/indices of columns to be retained.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing columns.

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
        aln.data.retain_cols(positions)
        aln.column_metadata = aln.column_metadata.iloc[positions]

        # # Add to history
        # add_to_history(
        #     self._instance, '.col.retain', positions,
        #     copy=copy,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def reorder(self, position_list, copy=False, **kwargs):
        """Reorders columns according the specified list of positions.
        
        Parameters
        ----------
        position_list : list of int
            Ordered list of position indices indicating the new order
            of columns.
        copy : bool, optional
            Whether to return a new copy of the reordered alignment, keeping the
            original intact, or reorder the alignment inplace. (default is
            False, reordering is done inplace)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after reordering
            columns.
        """
        # Check input
        if isinstance(position_list, list) and \
            sum((isinstance(pos, int) for pos in position_list)):
            if len(list) != self._instance.ncols:
                raise TypeError('length of position list must be equal to the '
                    'number of columns in the alignment: {} != {}'.format(
                        len(list), self._instance.ncols
                    ))
        else:
            raise TypeError('position list must be a list of int')

        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        aln.data.reorder_cols(position_list)
        aln.column_metadata = aln.column_metadata.iloc[position_list]

        # # Add to history
        # add_to_history(
        #     aln, '.col.reorder', positions,
        #     copy=copy,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def filter(self, function, copy=False, dry_run=False, inverse=False,
               chunk_size=1, **kwargs):
        """Returns the list of column positions where the given function
        is True.
        
        Parameters
        ----------
        function : callable
            Function used to evaluate each column. The function should expect a list of list of str as input and return a bool as output.
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
        chunk_size : int, optional
            Number of characters to group as one column. (default is 1)

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
        # Checks if chunk_size value is valid
        if chunk_size < 1:
            raise ValueError('chunk_size value must be greater than zero')

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

        # Get the list of positions based on the result of the
        # filtering function. Only columns that are True are recorded.
        # There are 2 different ways to create the positions list depending on
        # whether the chunk_size is 1 (default) or greater than 1.
        if chunk_size == 1:
            positions = [
            i for i, col in enumerate(aln.col.iter(chunk_size=chunk_size)) 
            if function(col)
        ]
        else:
            positions = list(itertools.chain(
                *[range(i*chunk_size, (i*chunk_size)+chunk_size)
                  for i in positions]))

        # Generate list of column ids that are "opposite" of what was given
        other_positions = aln.data.invert_cols(positions)

        # "dry_run" shows the columns that will are True or False based
        # on the given filtering function.
        # This prints out the number of columns in the True or False categories
        # and returns the lists of column ids that are classified as True or
        # False.
        if dry_run:
            parts = []
            parts.append('[{}]'.format(custom_title))
            parts.append('{} = {}/{}'.format(
                custom_class_true, len(positions), aln.ncols))
            parts.append('{} = {}/{}'.format(
                custom_class_false, len(other_positions), aln.ncols))
            print('\n'.join(parts))
            return {
                True: positions,
                False: other_positions
            }

        # By default, the filter method will keep columns that are
        # True according to the filter function, and will remove
        # columns that are False.
        # However, if `inverse` is True, the filter method will do the
        # opposite. It will keep columns that are False and will remove
        # columsn that are True
        if inverse:
            aln.col.remove(positions, _record_history=False)
        else:
            aln.col.retain(positions, _record_history=False)
            
        # # Add to history
        # func_sig = function.__qualname__ + \
        #     repr(inspect.signature(function)) \
        #         .lstrip('<Signature ').rstrip('>')
        # add_to_history(
        #     aln, '.col.filter', func_sig,
        #     copy=copy,
        #     dry_run=dry_run,
        #     inverse=inverse,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def has(self, query, case_sensitive=False, mode='any',
            step=None, chunk_size=None, **kwargs):
        """Return the list of positions that matches the given conditions.
        
        Parameters
        ----------
        query : str
            Query string
        case_sensitive : bool, optional
            Whether or not to look for a case-sensitive match. (default is False, matches are not case sensitive)
        mode : str, optional
            Specifies whether "any" or "all"-type matching is performed.
            In "any" mode, the column is a match if any of the items matches
            the query. In "all" mode, the column is a match if and only if
            all items match the query. (default is "any")
        step : int, optional
            Number of characters to skip. (default is 1)
        chunk_size : int, optional
            Number of characters to group as one column. (default is 1)

        Returns
        -------
        list
            List of column positions that matches the query.
        """
        # Check input
        if step and chunk_size:
            raise ValueError(
                'skip_n and chunk_size cannot be used simultaneously')
        if step is None:
            if chunk_size is None:
                skip_n = 1
            else:
                skip_n = chunk_size
        if chunk_size is None:
            chunk_size = 1
        
        # Pass inputs to Rust-backed function and output
        return self._instance.data.has(
            query, case_sensitive, mode, step, chunk_size)

    def map(self, function, step=None, chunk_size=None):
        """Maps a function to the sequence matrix column-wise.
        
        Parameters
        ----------
        function : callable
            Function to be mapped. The function should expect a list of str as input.
        step : int, optional
            Number of characters to skip. (default is None)
        chunk_size : int, optional
            Number of characters to group as one column. (default is None)

        Yields
        -------
        object
            Each column is evaluated using the given function

        """
        for col in self.iter(step=step, chunk_size=chunk_size, lazy=True):
            yield function(col)

    def iter(self, step=None, chunk_size=None, lazy=False):
        """Iterates over the sequence matrix column-wise.
        Returns a list of list of str, the inner list representing a column.

        Parameters
        ----------
        step : int, optional
            Number of characters to skip. (default is None)
        chunk_size : int, optional
            Number of characters to group as one column. (default is None)
        lazy : bool, optional
            If True, uses lazy execution (saves memory), otherwise uses eager execution. (default is False, uses eager execution)

        Yields
        -------
        list of str
            Each column is represented as a list of str whose order is based
            on the ordering of samples in the alignment.

        """
        # Check input type
        if step is not None:
            if not isinstance(step, int):
                raise ValueError(
                    '`step` must be None or int: {} ({})'.format(
                        step, type(step)
                    ))
        if chunk_size is not None:
            if not isinstance(chunk_size, int):
                raise ValueError(
                    '`chunk_size` must be None or int: {} ({})'.format(
                        chunk_size, type(chunk_size)
                    ))
        # Check values
        if step < 1:
            raise ValueError('`step` must be greater than zero: {}'.format(
                step
            ))
        if chunk_size < 1:
            raise ValueError(
                '`chunk_size` must be greater than zero: {}'.format(
                    chunk_size
                ))
        if step > 1 and chunk_size > 1:
            raise ValueError(
                '`step` and `chunk_size` cannot be used simultaneously')
        
        # Initialize values
        # Default for step and chunk_size is None
        # The values are changed depending on how step and chunk_size are
        # set by the user.
        # Note that both step and chunk_size cannot be set simultaneously by
        # the user.

        # step is not set, chunk_size is set
        if step is None:
            if chunk_size is None:
                # If both step and chunk_size are None, then step is 1
                step = 1
            else:
                # If step is None but chunk_size is specified, step
                # adopts the value of chunk_size to get consecutive
                # columns.
                step = chunk_size
        # chunk_size is not set
        if chunk_size is None:
            chunk_size = 1
        cnt = 0
        col_range = range(0, self._instance.ncols - (chunk_size-1), step)

        # iter method offers two ways to iterate: lazy and eager
        # In lazy execution, the function uses yield to return a copy of
        # the current column. Either get_chunk or get_col is used depending
        # on whether the chunk_size is specified or not
        if lazy:
            for i in col_range:
                if chunk_size == 1:
                    yield self._instance.data.get_col(i)
                else:
                    yield self._instance.data.get_chunk(i, chunk_size)
        # In eager execution, the function transforms the sequence matrix
        # into a list of list of str column-wise using get_chunks or get_cols
        # depending on whether the chunk_size is specified or not.
        # Then the list of list of str is iterated, using yield to return
        # each column (list of str) one by one.
        else:
            indices = list(col_range)
            if chunk_size == 1:
                for col in self._instance.data.get_cols(indices):
                    yield col
            else:
                for col in self._instance.data.get_chunks(indices, chunk_size):
                    yield col

    def reset_index(self, copy=False, drop=False, **kwargs):
        """Resets the column index to the default integer index.
        
        Parameters
        ----------
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)
        drop : bool, optional
            If True, do not try to insert the original index into dataframe
            columns. (default is False, the original index is inserted as a
            column named `index`)
        
        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing columns that evaluated False. Note that this returns the whole
            Alignment object and not only the pandas DataFrame containing
            column metadata.

        """
        aln = self._instance
        if copy is True:
            aln = self._instance.copy()
        aln.column_metadata.reset_index(drop=drop, inplace=True)

        # # Add to history
        # add_to_history(
        #     aln, '.col.reset_index',
        #     copy=copy,
        #     **kwargs
        # )

        if copy is True:
            return aln

    def add_metadata(self, metadata, name=None, copy=False, **kwargs):
        """Adds a new category to the column metadata. This adds a column
        to the column metadata DataFrame.
        
        Parameters
        ----------
        metadata : list, dict, pandas.Series or pandas.DataFrame
            Metadata to be added.
        name : str, optional
            Name of the new metadata category.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)

        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after adding new column metadata categories.
        
        """
        raise NotImplementedError()

    def remove_metadata(self, name, copy=False, **kwargs):
        """Removes one or more categories from the column metadata. This removes
        columns from the column metadata DataFrame.
        
        Parameters
        ----------
        name : str or list of str
            Name/s of the new metadata categories.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)

        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after removing column metadata categories.
        
        """
        raise NotImplementedError()

    def replace_metadata(self, name, metadata, copy=False):
        """Replaces metadata in the given column metadata category.
        
        Parameters
        ----------
        name : str
            Name of the metadata category. This is also the column name
            in the underlying column metadata DataFrame.
        metadata: list
            List of metadata to replace existing information.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace. (default is False,
            editing is done inplace)

        Returns
        -------
        Alignment
            When `copy` is True, returns the edited alignment after replacing
            the metadata in the specified column metadata category.
        
        """
        raise NotImplementedError()

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
