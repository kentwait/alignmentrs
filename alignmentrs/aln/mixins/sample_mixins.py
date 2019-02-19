from .functions import *


__all__ = ['SamplePropsMixin', 'SampleAlnMixin']

class SamplePropsMixin:
    @property
    def nsamples(self):
        """int: Returns the number of samples in the alignment."""
        if not self.samples:
            return 0
        return self.samples.nrows

    @property
    def sample_ids(self):
        """list of str: Returns the list of sample identifiers."""
        if not self.samples:
            return []
        return self.samples.ids

    @property
    def sample_descriptions(self):
        """list of str: Returns the list of sample descriptions."""
        if not self.samples:
            return []
        return self.samples.descriptions

    @property
    def sample_sequences(self):
        """list of str: Returns the list of sample sequences."""
        if not self.samples:
            return []
        return self.samples.sequences


class MarkerPropsMixin:
    @property
    def nmarkers(self):
        """int: Returns the number of markers in the alignment."""
        if not self.markers:
            return 0
        return self.markers.nrows

    @property
    def marker_ids(self):
        """list of str: Returns the list of markers identifiers."""
        if not self.markers:
            return []
        return self.markers.ids

    @property
    def marker_descriptions(self):
        """list of str: Returns the list of marker descriptions."""
        if not self.markers:
            return []
        return self.markers.descriptions

    @property
    def marker_sequences(self):
        """list of str: Returns the list of marker sequences."""
        if not self.markers:
            return []
        return self.markers.sequences

class SampleAlnMixin:
    # Insert Methods
    # ==========================================================================

    # Insert/append one sample
    # ------------------------------
    def insert_sample(self, row, name, description, sequence, copy=False):
        """Inserts a sample entry at the specified position.

        Parameters
        ----------
        row : int
            Row position to insert the sample.
        name : str
            Sample identifier.
        description : str
            Sample description.
        sequence : str
            Sample sequence.
        copy : bool, optional
            Returns a new copy instead of inserting inplace.
            (default is False, operation is done inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        insert_items_from_lists(
            aln.samples, row, [name], [description], [sequence])
        if copy:
            return aln

    def append_sample(self, name, description, sequence, copy=False):
        """Appends a sample entry after the last sample in the alignment.

        Parameters
        ----------
        name : str
            Sample identifier.
        description : str
            Sample description.
        sequence : str
            Sample sequence.
        copy : bool, optional
            Returns a new copy instead of inserting inplace.
            (default is False, operation is done inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        append_items_from_lists(
            aln.samples, [name], [description], [sequence])
        if copy:
            return aln


    # Insert/append many samples
    # ------------------------------
    # TODO: Use record for insert/append many
    def insert_samples(self, row, records, copy=False):
        raise NotImplementedError()

    def append_samples(self, records, copy=False):
        raise NotImplementedError()

    def insert_samples_from_lists(self, row, ids, descriptions, sequences, 
                                  copy=False):
        """Inserts sample entries at the specified position.

        Parameters
        ----------
        row : int
            Row position to insert the samples.
        ids : list of str
            List of sample identifiers to insert. The order should correspond
            to the order of ther other lists.
        descriptions : list of str
            List of sample descriptions to insert. The order should correspond
            to the order of ther other lists.
        sequences : list of str
            List of sample sequneces to insert. The order should correspond
            to the order of ther other lists.
        copy : bool, optional
            Returns a new copy instead of inserting inplace.
            (default is False, operation is done inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        insert_items_from_lists(aln.samples, row, ids, descriptions, sequences)
        if copy:
            return aln

    def append_samples_from_lists(self, ids, descriptions, sequences, copy=False):
        """Appends sample entries after the last sample in the alignment.

        Parameters
        ----------
        ids : list of str
            List of sample identifiers to insert. The order should correspond
            to the order of ther other lists.
        descriptions : list of str
            List of sample descriptions to insert. The order should correspond
            to the order of ther other lists.
        sequences : list of str
            List of sample sequneces to insert. The order should correspond
            to the order of ther other lists.
        copy : bool, optional
            Returns a new copy instead of inserting inplace.
            (default is False, operation is done inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        append_items_from_lists(aln.samples, ids, descriptions, sequences)
        aln.samples.append_rows(ids, descriptions, sequences)
        if copy:
            return aln


    # Setter methods
    # ==========================================================================
    def replace_samples(self, rows, sequences, copy=False):
        """Replaces sequences of the specified samples.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the samples to be replaced.
            If a list, length must match the number of given sequences.
        sequences : str
            A single sequence or a list of sequences.
            If `rows` is an int/str, `sequences` must be a str.
            If `rows` is a list, `sequences` must list and the number of sequences
            must match the length of `rows`.
        copy : bool, optional
            Returns a new copy instead of performing the
            replacement inplace. (default is False, operation is done
            inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        replace_items(aln.samples, rows, sequences)
        if copy:
            return aln

    def reorder_samples(self, rows, copy=False):
        """Reorders samples based on a list of identifiers or indices.

        Parameters
        ----------
        rows : list of int or list of str
            Order of the list specifies the new ordering of the samples.
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
        reorder_items(aln.samples, rows)
        if copy:
            return aln


    # Sample deleters
    # ------------------------------
    def remove_samples(self, rows, match_prefix=False, match_suffix=False,
                       copy=False):
        """Removes samples based a list of identifiers or indices.
        This is the opposite of the `retain_samples` method.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the samples to be removed.
        match_prefix : bool, optional
            Whether to interpret `rows` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `rows` as a suffix to match against
            the list of sample names. This parameter is considered
            only if match_prefix is False. (default is False)
        copy : bool, optional
            Returns a new copy instead of removing samples inplace.
            (default is False, operation is done inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        remove_items(aln.samples, rows, match_prefix=match_prefix,
                     match_suffix=match_suffix)
        if copy:
            return aln

    def retain_samples(self, rows, match_prefix=False, match_suffix=False, copy=False):
        """Removes samples based a list of identifiers or indices.
        This is the opposite of the `remove_samples` method.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the samples to be retained.
        match_prefix : bool, optional
            Whether to interpret `rows` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `rows` as a suffix to match against
            the list of sample names. This parameter is considered
            only if match_prefix is False. (default is False)
        copy : bool, optional
            Returns a new copy instead of performing the
            operation inplace. (default is False, operation is done
            inplace)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.copy() if copy else self
        retain_items(aln.samples, rows, match_prefix=match_prefix, 
                     match_suffix=match_suffix)
        if copy:
            return aln


    # Iterators
    # ==========================================================================
    def iter_sample_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the sample alignment. Excludes markers.

        Parameters
        ----------
        start : int, optional
            Starting position. (default is 0)
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
        for v in iter_aln_sites(self.samples):
            yield v

    def iter_samples(self):
        """Iterates over samples in the alignment.
        Excludes markers.

        Yields
        ------
        str

        """
        for v in iter_aln(self.samples):
            yield v

    def iter_sample_records(self):
        """Iterates over samples in the alignment, returning a Record object.
        Excludes markers.

        Yields
        ------
        Record

        """
        for v in iter_aln_records(self.samples):
            yield v

    def get_subset(self, rows=None, cols=None):
        """Returns a subset of the alignment based on the given set of
        samples, markers and sites.

        Parameters
        ----------
        rows : int, list of int, or None
            An int/str/list specifying the samples to be included.
            If None, all samples will be included in the subset.
        cols : int, list of int, or None
            int, or list specifying the sites to be included.
            If None, all sites will be included in the subset.

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.
        ValueError
            marker_ids is specified but the alignment has no
            marker sequences.

        Returns
        -------
        Alignment
            New alignment object containing the subset of sample and
            markers rows, and site columns.
            This subset is a deep copy of the original alignment and
            will not be affect by changes made in the original.

        """
        if rows is None:
            return subset(self, cols)
        return subset(self, cols, samples=rows)
