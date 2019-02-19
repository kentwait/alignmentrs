from copy import deepcopy
import itertools

from libalignmentrs.record import Record
from libalignmentrs.position import BlockSpace
from .functions import subset


__all__ = ['PropsMixin', 'AlnMixin']

class PropsMixin:
    @property
    def nrows(self):
        """int: Returns the number of rows in the alignment."""
        return sum(
            (self.__getattribute__(m).nrows for m in self.__class__.members))

    @property
    def ncols(self):
        """int: Returns the number of sites in the alignment."""
        return self.samples.ncols

    @property
    def ids(self):
        """list of str: Returns the list of identifiers."""
        generator = (
            self.__getattribute__(m).ids for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def descriptions(self):
        """list of str: Returns the list of descriptions."""
        generator = (
            self.__getattribute__(m).descriptions 
            for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def sequences(self):
        """list of str: Returns the list of sequences."""
        generator = (
            self.__getattribute__(m).sequences 
            for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def coordinates(self):
        """list of int: Returns the list of coordinates
        associated to each column in the alignment."""
        return self._linspace.to_arrays()[0]


class AlnMixin:
    # Getters
    # ------------------------------
    def get_rows(self, rows, match_prefix=False, match_suffix=False):
        """Returns a list of records containing only the specified samples.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the samples to be retrieved.
        match_prefix : bool, optional
            Whether to interpret `rows` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `rows` as a suffix to match against
            the list of sample names. This parameter is considered
            only if match_prefix is False. (default is False)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        List of Record
            Returns matches as a list of Record object.

        """
        records = []
        # Offsets the int indices
        ro = 0
        for member in self.__class__.members:
            if isinstance(rows, int):
                rows = [row-ro for row in rows if row-ro > 0]
                records += self.__getattribute__(member).get_rows(rows)
            elif isinstance(rows, str):
                if match_prefix:
                    records += self.__getattribute__(member) \
                        .get_rows_by_prefix([rows])
                elif match_suffix:
                    records +=  self.__getattribute__(member) \
                        .get_rows_by_suffix([rows])
                else:
                    records +=  self.__getattribute__(member) \
                        .get_rows_by_name([rows])
            elif isinstance(rows, list) and \
                sum((isinstance(j, int) for j in rows)):
                rows = [row-ro for row in rows if row-ro > 0]
                records +=  self.__getattribute__(member).get_rows(rows)
            elif isinstance(rows, list) and \
                sum((isinstance(j, str) for j in rows)):
                if match_prefix:
                    records +=  self.__getattribute__(member) \
                        .get_rows_by_prefix(rows)
                elif match_suffix:
                    records +=  self.__getattribute__(member) \
                        .get_rows_by_suffix(rows)
                else:
                    records +=  self.__getattribute__(member) \
                        .get_rows_by_name(rows)
            else:
                raise TypeError('rows must be an int, str, list of int, or list of str.')
            ro += self.__getattribute__(member).nrows()
        return records

    def get_cols(self, cols):
        """Returns a list of records representing each alignment column specified by index/list of indices.

        Parameters
        ----------
        cols : int or list of int
            An int/list specifying the sites to be retrieved.

        Returns
        -------
        Alignment
            Creates a new Alignment object containing the specified the
            specified sites.
            This subset is a deep copy of the original alignment
            and will be independent of changes made in the original.

        """
        member_records = []
        # Offsets the int indices
        ro = 0
        for member in self.__class__.members:
            records = []
            if isinstance(cols, int):
                records += self.__getattribute__(member).get_cols([cols])
            elif isinstance(cols, list) and \
                sum((isinstance(j, int) for j in cols)):
                cols = [row-ro for row in cols if row-ro > 0]
                records +=  self.__getattribute__(member).get_cols(cols)
            else:
                raise TypeError('cols must be an int, or list of int.')
            member_records.append(records)
        # Combine results
        records = []
        for i in range(len(cols)):
            rid = member_records[0][i].id
            description = member_records[0][i].description
            sequence = ''.join([member_records[m][i].sequence
                                for m in range(len(member_records))])
            record = Record(rid, description, sequence)
            records.append(record)
        return records

    def reset_coordinates(self, start=0, stop=None, state=1):
        """Resets the coordinates of the alignment columns.  

        Parameters
        ----------
        start : int, optional
            Starting coordinate value (default is 0)
        stop : int, optional
            End coordinate value.
            This value is not part of linear space representing the alignment.
            (default is None, ending value becomes the sum of start and
            the number of sites)
        state : int, optional
            State annotation. (default is 1)

        """
        start = start if start is not None else 0
        stop = stop if stop is not None else start + self.samples.ncols
        state = state if state is not None else 1
        self._linspace: BlockSpace = BlockSpace(start, stop, state)


    # Deleters
    # ------------------------------
    def remove_cols(self, cols, copy=False):
        """Removes sites based on a list of column numbers.
        This is the opposite of the `retain_cols` method.

        Parameters
        ----------
        cols : int, or list of int
            An int/list specifying the sites to be removed.
        copy : bool, optional
            Returns a new copy instead of removing sites inplace.
            (default is False, operation is done inplace)

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        # Check type of cols, and convert if necessary
        if isinstance(cols, int):
            cols = [cols]
        # Perform removal inplace
        columns = None
        for member in self.__class__.members:
            aln.__getattribute__(member).remove_cols(cols)
            if columns is None:
                columns = aln.__getattribute__(member).ncols
            else:
                if columns != aln.__getattribute__(member).ncols:
                    raise ValueError('Column number mismatch')
        aln._linspace.remove(cols)
        if copy:
            return aln

    def retain_cols(self, cols, copy=False):
        """Keeps sites based on a list of column numbers.
        This is the opposite of the `remove_cols` method.

        Parameters
        ----------
        cols : int or list of int
            An int/list specifying the sites to be retained.
        copy : bool, optional
            Returns a new copy instead of performing the operation inplace.
            (default is False, operation is done inplace)

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        # Check type of cols, and convert if necessary
        if isinstance(cols, int):
            cols = [cols]
        # Perform removal inplace
        columns = None
        for member in self.__class__.members:
            aln.__getattribute__(member).retain_cols(cols)
            if columns is None:
                columns = aln.__getattribute__(member).ncols
            else:
                if columns != aln.__getattribute__(member).ncols:
                    raise ValueError('Column number mismatch')
        aln._linspace.retain(cols)
        if copy:
            return aln


    # Iterators
    # ==========================================================================
    def iter_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the alignment.

        Parameters
        ----------
        start : int, optional
            Starting column position. (default is 0)
        stop : [type], optional
            Stopping column position.
            If None, the iterator will continue until the end of the alignment.
            (default is None)
        size : int, optional
            Number of characters to yield at each iteration.
            For single characters, `size` = 1.
            For codons, `size` = 3. (default is 1)

        Raises
        ------
        ValueError
            If the alignment cannot be cleanly cut up into the specified
            chunk size (`ncols` not divisible be `size`),
            a ValueError is raised.

        Yields
        ------
        list of str
            List of sequences representing a site or chunk of the alignment.

        """
        if stop is None:
            stop = self.ncols
        if (stop - start) % size != 0:
            raise ValueError('Alignment cannot be completely divided into '
                             'chucks of size {}'.format(size))
        for i in range(start, stop, size):
            results = []
            for member in self.__class__.members:
                results += [s[i:i+size] for s in 
                            self.__getattribute__(member).sequences]
                yield results

    def iter_rows(self):
        """Iterates over all rows in the alignment.

        Yields
        ------
        str

        """
        for member in self.__class__.members:
            for sequence in self.__getattribute__(member).sequences:
                yield sequence
    
    def iter_records(self):
        """Iterates over all rows in the alignment, returning a Record object.

        Yields
        ------
        Record

        """
        for member in self.__class__.members:
            for vals in zip(self.__getattribute__(member).ids,
                            self.__getattribute__(member).descriptions,
                            self.__getattribute__(member).sequences):
                yield Record(*vals)
