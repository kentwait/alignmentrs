from .functions import *


__all__ = ['MarkerPropsMixin', 'MarkerAlnMixin']

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


class MarkerAlnMixin:
    # Insert Methods
    # ==========================================================================

    # Insert/append one marker
    # ------------------------------
    def insert_marker(self, row, name, description, sequence, copy=False):
        """Inserts a marker entry at the specified position.

        Parameters
        ----------
        row : int
            Row position to insert the marker.
        name : str
            marker identifier.
        description : str
            marker description.
        sequence : str
            marker sequence.
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
            aln.markers, row, [name], [description], [sequence])
        if copy:
            return aln

    def append_marker(self, name, description, sequence, copy=False):
        """Appends a marker entry after the last marker in the alignment.

        Parameters
        ----------
        name : str
            marker identifier.
        description : str
            marker description.
        sequence : str
            marker sequence.
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
            aln.markers, [name], [description], [sequence])
        if copy:
            return aln


    # Insert/append many markers
    # ------------------------------
    # TODO: Use record for insert/append many
    def insert_markers(self, row, records, copy=False):
        raise NotImplementedError()

    def append_markers(self, records, copy=False):
        raise NotImplementedError()

    def insert_markers_from_lists(self, row, ids, descriptions, sequences, 
                                  copy=False):
        """Inserts marker entries at the specified position.

        Parameters
        ----------
        row : int
            Row position to insert the markers.
        ids : list of str
            List of marker identifiers to insert. The order should correspond
            to the order of ther other lists.
        descriptions : list of str
            List of marker descriptions to insert. The order should correspond
            to the order of ther other lists.
        sequences : list of str
            List of marker sequneces to insert. The order should correspond
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
        insert_items_from_lists(aln.markers, row, ids, descriptions, sequences)
        if copy:
            return aln

    def append_markers_from_lists(self, ids, descriptions, sequences,
                                  copy=False):
        """Appends marker entries after the last marker in the alignment.

        Parameters
        ----------
        ids : list of str
            List of marker identifiers to insert. The order should correspond
            to the order of ther other lists.
        descriptions : list of str
            List of marker descriptions to insert. The order should correspond
            to the order of ther other lists.
        sequences : list of str
            List of marker sequneces to insert. The order should correspond
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
        append_items_from_lists(aln.markers, ids, descriptions, sequences)
        if copy:
            return aln


    # Setter methods
    # ==========================================================================
    def replace_markers(self, rows, sequences, copy=False):
        """Replaces sequences of the specified markers.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the markers to be replaced.
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
        replace_items(aln.markers, rows, sequences)
        if copy:
            return aln

    def reorder_markers(self, rows, copy=False):
        """Reorders markers based on a list of identifiers or indices.

        Parameters
        ----------
        rows : list of int or list of str
            Order of the list specifies the new ordering of the markers.
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
        reorder_items(aln.markers, rows)
        if copy:
            return aln


    # marker deleters
    # ------------------------------
    def remove_markers(self, rows, match_prefix=False, match_suffix=False,
                       copy=False):
        """Removes markers based a list of identifiers or indices.
        This is the opposite of the `retain_markers` method.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the markers to be removed.
        match_prefix : bool, optional
            Whether to interpret `rows` as a prefix to match against
            the list of marker names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `rows` as a suffix to match against
            the list of marker names. This parameter is considered
            only if match_prefix is False. (default is False)
        copy : bool, optional
            Returns a new copy instead of removing markers inplace.
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
        remove_items(aln.markers, rows, match_prefix=match_prefix, 
                     match_suffix=match_suffix)
        if copy:
            return aln

    def retain_markers(self, rows, match_prefix=False, match_suffix=False,
                       copy=False):
        """Removes markers based a list of identifiers or indices.
        This is the opposite of the `remove_markers` method.

        Parameters
        ----------
        rows : int, str, list of int, or list of str
            An int/str/list specifying the markers to be retained.
        match_prefix : bool, optional
            Whether to interpret `rows` as a prefix to match against
            the list of marker names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `rows` as a suffix to match against
            the list of marker names. This parameter is considered
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
        retain_items(aln.markers, rows, match_prefix=match_prefix, 
                     match_suffix=match_suffix)
        if copy:
            return aln


    # Iterators
    # ==========================================================================
    def iter_marker_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the marker alignment. Excludes markers.

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
        for v in iter_aln_sites(self.markers):
            yield v

    def iter_markers(self):
        """Iterates over markers in the alignment.
        Excludes markers.

        Yields
        ------
        str

        """
        for v in iter_aln(self.markers):
            yield v

    def iter_marker_records(self):
        """Iterates over markers in the alignment, returning a Record object.
        Excludes markers.

        Yields
        ------
        Record

        """
        for v in iter_aln_records(self.markers):
            yield v

    def get_subset(self, rows, cols):
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
        return subset(self, cols, markers=rows)
