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
            Column index/indices of column/s to return.
        
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
            Column index/indices of columns to be removed.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace.
        
        Returns
        -------
        Alignment
            Returns the edited alignment after removing the specified
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
            Column index/indices of columns to be retained.
        copy : bool, optional
            Whether to return a new copy of the edited alignment, keeping the
            original intact, or edit the alignment inplace.
        
        Returns
        -------
        Alignment
            Returns the edited alignment after removing columns.

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
            Ordered list of column indices indicating the new order of columns.
        copy : bool, optional
            Whether to return a new copy of the reordered alignment, keeping the
            original intact, or reorder the alignment inplace.
        
        Raises
        ------
        TypeError
            [description]
        
        Returns
        -------
        [type]
            [description]
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
            List of column positions that match the query.
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

    def drop(self, value, case_sensitive=False, copy=False, dry_run=False,
             mode='any', **kwargs):
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

    def map(self, function, skip_n=None, chunk_size=None, lazy=False):
        for col in self.iter(skip_n=skip_n, chunk_size=chunk_size, lazy=lazy):
            yield function(col)

    def iter(self, skip_n=None, chunk_size=None, lazy=False):
        cnt = 0
        if skip_n and chunk_size:
            raise ValueError(
                'skip_n and chunk_size cannot be used simultaneously')
        if skip_n is None:
            if chunk_size is None:
                skip_n = 1
            else:
                skip_n = chunk_size
        if chunk_size is None:
            chunk_size = 1

        col_range = range(0, self._instance.ncols-(chunk_size-1), skip_n)
        if lazy:
            for i in col_range:
                if chunk_size == 1:
                    yield self._instance.data.get_col(i)
                else:
                    yield self._instance.data.get_chunk(i, chunk_size)
        else:
            indices = list(col_range)
            if chunk_size == 1:
                for col in self._instance.data.get_cols(indices):
                    yield col
            else:
                for col in self._instance.data.get_chunks(indices, chunk_size):
                    yield col

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
