import numpy as np
from libalignmentrs.alignment import BaseAlignment, fasta_file_to_basealignments
from libalignmentrs.position import BlockSpace
from alignmentrs.util import parse_comment_list


__all__ = ['Alignment']


class Alignment:
    """Represents a multiple sequence alignment.

    Alignment encapsulates information included in the FASTA format:
    sequence names/ids, descriptions, and sequences.

    Additionally, Alignment can track site order using
    block metadata.

    Attributes
    ----------
    name : str
        Name of the alignment.
    samples : BaseAlignment
        Alignment of sample sequences
    markers : BaseAlignment.
        Alignment of non-sample sequences. This is metadata
        stored as a row in the alignment that describes some
        kind of site-specific information.
    blocklists : list of Block
        List that stores positional information of the alignment since
        tracking was initiated.

    """

    def __init__(self, name, sample_alignment, marker_alignment, **kwargs):
        """Creates a new Alignment object from sample and marker alignments.

        Parameters
        ----------
        name : str
            Name of the alignment.
        sample_alignment : BaseAlignment
            Alignment of sample sequences.
        marker_alignment: BaseAlignment
            Alignment of non-sample sequences. These sequences are
            usually metadata to indicate site-specific information.
            If None, an empty BaseAlignment replaces the None value.

        Raises
        ------
        ValueError
            Alignment is instantiated with an empty sample alignment, or
            instantiated with sample and marker alingments of unequal number of
            sites.

        """
        if not sample_alignment:
            raise ValueError(
                'Cannot create an Alignment using an empty '
                'sample_alignment BaseAlignment.')
        if marker_alignment and \
           sample_alignment.nsites != marker_alignment.nsites:
            raise ValueError(
                'Number of sites in sample and marker alignments are '
                'not equal: {} != {}'.format(
                    sample_alignment.nsites, marker_alignment.nsites))
        self.name = name
        self.samples: BaseAlignment = sample_alignment
        self.markers: BaseAlignment = marker_alignment \
            if marker_alignment else BaseAlignment([], [], [])
        if '_linspace' in kwargs.keys():
            self._linspace: BlockSpace = kwargs['_linspace']
        else:
            start = kwargs['linspace_default_start'] \
                    if 'linspace_default_start' in kwargs.keys() else 0
            stop = kwargs['linspace_default_stop'] \
                    if 'linspace_default_stop' in kwargs.keys() else \
                    start + self.samples.nsites
            state = kwargs['linspace_default_state'] \
                    if 'linspace_default_state' in kwargs.keys() else \
                    1
            self._linspace: BlockSpace = BlockSpace(start, stop, state)

    # Properties
    # ==========================================================================

    # General properties
    # ------------------------------
    @property
    def nrows(self):
        """int: Returns the number of rows in the alignment."""
        nmarkers = self.markers.nrows if self.nmarkers else 0
        return self.samples.nrows + nmarkers

    @property
    def nsites(self):
        """int: Returns the number of sites in the alignment."""
        return self.samples.nsites

    @property
    def coordinates(self):
        """list of int: Returns the coordinates of the current alignment
        relative to positions when the alignment was created/loaded."""
        return self._linspace.to_arrays()[0]

    # Sample properties
    # ------------------------------
    @property
    def nsamples(self):
        """int: Returns the number of samples in the alignment."""
        return self.samples.nrows

    @property
    def sample_ids(self):
        """list of str: Returns the list of sample sequences."""
        return self.samples.ids

    @property
    def sample_descriptions(self):
        """list of str: Returns the list of sample sequences."""
        return self.samples.descriptions

    @property
    def sample_sequences(self):
        """list of str: Returns the list of sample sequences."""
        return self.samples.sequences

    # Marker properties
    # ------------------------------
    @property
    def nmarkers(self):
        """int: Returns the number of markers in the alignment."""
        if not self.markers:
            return 0
        return self.markers.nrows

    @property
    def marker_ids(self):
        """list of str: Returns the list of sample sequences."""
        if not self.markers:
            return []
        return self.markers.ids

    @property
    def marker_descriptions(self):
        """list of str: Returns the list of sample sequences."""
        if not self.markers:
            return []
        return self.markers.descriptions

    @property
    def marker_sequences(self):
        """list of str: Returns the list of sample sequences."""
        if not self.markers:
            return []
        return self.markers.sequences


    # Getter methods
    # ==========================================================================

    # General getters
    # ------------------------------
    @classmethod
    def subset(cls, aln, sample_ids=None, marker_ids=None, sites=None):
        """Returns a subset of the alignment by samples, markers and sites.

        Parameters
        ----------
        aln : Alignment
        sample_ids : int, list of int, or None
            int, str or list specifying the samples to be included.
            If None, all samples will be included in the subset.
        marker_ids : int, list of int, or None
            int, str or list specifying the markers to be included.
            Row indices for markers in the alignment.
            If None, all markers will be included in the subset.
        sites : int, list of int, or None
            int, or list specifying the sites to be included.
            If None, all sites will be included in the subset.

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.
        ValueError
            marker_ids is specified by aln.markers is empty.

        Returns
        -------
        Alignment
            New alignment object containing the subset of sample and
            markers rows, and site columns.
            This subset is a deep copy of the original alignment and
            will not be affect by changes made in the original.

        """
        # Checks the value of sample_ids and converts if necessary.
        if sample_ids is None:
            sample_ids = list(range(0, aln.nsamples))
        elif isinstance(sample_ids, int):
            sample_ids = [sample_ids]
        elif isinstance(sample_ids, str):
            sample_ids = aln.samples.row_names_to_ids([sample_ids])
        elif (isinstance(sample_ids, list) and
              sum((isinstance(j, int) for j in sample_ids))):
            pass
        elif (isinstance(sample_ids, list) and
              sum((isinstance(j, str) for j in sample_ids))):
            sample_ids = aln.samples.row_names_to_ids(sample_ids)
        else:
            raise TypeError('sample_ids must be an int, str, list of int, '
                            'or list of str.')
        # Check if marker_ids is not None and checks if markers exist
        if marker_ids and not aln.markers:
            raise ValueError('Markers are not present in this alignment.')
        # Checks the value of marker_ids and converts if necessary.
        if marker_ids is None:
            marker_ids = list(range(0, aln.nsites))
        elif isinstance(marker_ids, int):
            marker_ids = [marker_ids]
        elif isinstance(marker_ids, str):
            marker_ids = aln.samples.row_names_to_ids([marker_ids])
        elif (isinstance(marker_ids, list) and
              sum((isinstance(j, int) for j in marker_ids))):
            pass
        elif (isinstance(marker_ids, list) and
              sum((isinstance(j, str) for j in marker_ids))):
            marker_ids = aln.samples.row_names_to_ids(marker_ids)
        else:
            raise TypeError('marker_ids must be an int, str, list of int, '
                            'or list of str.')
        # Checks the value of sites and converts if necessary.
        if sites is None:
            sites = list(range(0, aln.nsites))
        elif isinstance(sites, int):
            sites = [sites]
        elif (isinstance(sites, list) and
              sum((isinstance(j, int) for j in sites))):
            pass
        else:
            raise TypeError('Sites must be an int, or list of int.')
        # Create new BaseAlignments for sample and marker,
        # if it exists in the original
        sample_aln = aln.samples.subset(sample_ids, sites)
        marker_aln = aln.samples.subset(marker_ids, sites) if aln.markers else \
                     None
        return cls(aln.name, sample_aln, marker_aln,
                   _linspace=aln._linspace.extract(sites))

    def get_sites(self, i):
        """Returns a new alignment containing only the sites specified
        by the given list of column numbers.

        Parameters
        ----------
        i : int or list of int
            int or list specifying the sites to be retrieved.

        Returns
        -------
        Alignment
            Creates a new Alignment object containing the specified the
            specified sites.
            This subset is a deep copy of the original alignment
            and will be independent of changes made in the original.

        """
        return self.__class__.subset(self, sites=i)

    # Sample getters
    # ------------------------------
    def get_samples(self, i, match_prefix=False, match_suffix=False):
        """Returns a list of sequence strings containing only the samples
        specified by the index.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the samples to be retrieved.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
            the list of sample names. This parameter is considered
            only if match_prefix is False. (default is False)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        Alignment
            Creates a new Alignment object containing the specified sample/s
            and markers are not included. This subset is a deep copy of the
            original slignment and will not be affect by changes made in the
            original.

        """
        # Call get_sample/s method for sample BaseAlignment depending on the
        # type of i
        sample_aln = None
        if isinstance(i, int):
            sample_aln = self.samples.get_rows([i])
        elif isinstance(i, str):
            if match_prefix:
                sample_aln = self.samples.get_rows_by_prefix([i])
            elif match_suffix:
                sample_aln = self.samples.get_rows_by_suffix([i])
            else:
                sample_aln = self.samples.get_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            sample_aln = self.samples.get_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                sample_aln = self.samples.get_rows_by_prefix(i)
            elif match_suffix:
                sample_aln = self.samples.get_rows_by_suffix(i)
            else:
                sample_aln = self.samples.get_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if sample_aln is None:
            raise ValueError('Value of `sample_aln` cannot be None.')
        return self.__class__(self.name, sample_aln, None,
                              _linspace=self._linspace.copy())

    # Marker getters
    # ------------------------------
    def get_markers(self, i, match_prefix=False, match_suffix=False):
        """Returns a list of sequence strings containing only the markers
        specified by the index.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the markers to be retrieved.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of markers names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
            the list of markers names. This parameter is considered
            only if match_prefix is False. (default is False)

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.

        Returns
        -------
        BaseAlignment
            Creates a new marker BaseAlignment object containing the
            specified sample/s.
            This subset is a deep copy of the original marker BaseAlignment
            object and will not be affect by changes made in the original.

        """
        # Call get_sample/s method for sample BaseAlignment depending on the
        # type of i
        if isinstance(i, int):
            return self.markers.get_rows([i])
        elif isinstance(i, str):
            if match_prefix:
                return self.markers.get_rows_by_prefix([i])
            elif match_suffix:
                return self.markers.get_rows_by_suffix([i])
            else:
                return self.markers.get_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            return self.markers.get_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                return self.markers.get_rows_by_prefix(i)
            elif match_suffix:
                return self.markers.get_rows_by_suffix(i)
            else:
                return self.markers.get_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')


    # Insert Methods
    # ==========================================================================

    # Insert/append one sample
    # ------------------------------
    def insert_sample(self, i, name, description, sequence, copy=False):
        """Inserts new sequences in the alignment matrix at the specified
        row position inplace.

        Parameters
        ----------
        i : int
            Row position to insert the samples into the sample alignment.
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
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
        self.insert_samples_from_lists(
            i, [name], [description], [sequence], copy=copy
        )

    def append_sample(self, name, description, sequence, copy=False):
        """Inserts new sequences after the last row of the alignment matrix
        inplace. This increases the total number of samples.

        Parameters
        ----------
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
        copy : bool, optional
            Returns a new copy instead of appending the samples inplace.
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
        self.append_samples_from_lists(
            [name], [description], [sequence], copy=copy
        )

    # Insert/append many samples
    # ------------------------------
    def insert_samples_from_lists(self, i, ids, descriptions, samples, copy=False):
        """Inserts new sequences in the alignment matrix at the specified
        row position inplace.

        Parameters
        ----------
        i : int
            Row position to insert the samples into the sample alignment.
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if not(isinstance(ids, list) and
               sum((isinstance(j, str) for j in ids))):
            raise TypeError('ids must be a list of str.')
        if not(isinstance(descriptions, list) and
               sum((isinstance(j, str) for j in descriptions))):
            raise TypeError('descriptions must be a list of str.')
        if not(isinstance(samples, list) and
               sum((isinstance(j, str) for j in samples))):
            raise TypeError('samples must be a list of str.')
        aln.samples.insert_rows(i, ids, descriptions, samples)
        if copy:
            return aln

    def append_samples_from_lists(self, ids, descriptions, samples, copy=False):
        """Inserts new sequences after the last row of the alignment matrix
        inplace. This increases the total number of samples.

        Parameters
        ----------
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
        copy : bool, optional
            Returns a new copy instead of appending the samples inplace.
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if not(isinstance(ids, list) and
               sum((isinstance(j, str) for j in ids))):
            raise TypeError('ids must be a list of str.')
        if not(isinstance(descriptions, list) and
               sum((isinstance(j, str) for j in descriptions))):
            raise TypeError('descriptions must be a list of str.')
        if not(isinstance(samples, list) and
               sum((isinstance(j, str) for j in samples))):
            raise TypeError('samples must be a list of str.')
        aln.samples.append_rows(ids, descriptions, samples)
        if copy:
            return aln

    # Insert/append one marker
    # ------------------------------
    def insert_marker(self, i, name, description, sequence, copy=False):    
        """Inserts new sequences in the alignment matrix at the specified
        row position inplace.

        Parameters
        ----------
        i : int
            Row position to insert the samples into the sample alignment.
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
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
        self.insert_markers_from_lists(
            i, [name], [description], [sequence], copy=copy
        )

    def append_marker(self, name, description, sequence, copy=False):
        """Inserts new sequences after the last row of the alignment matrix
        inplace. This increases the total number of samples.

        Parameters
        ----------
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        samples : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
        copy : bool, optional
            Returns a new copy instead of appending the samples inplace.
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
        self.append_markers_from_lists(
            [name], [description], [sequence], copy=copy
        )

    # Insert/append many markers
    # ------------------------------
    def insert_markers_from_lists(self, i, ids, descriptions, markers, copy=False):
        """Inserts new sequences in the alignment matrix at the specified
        row position inplace.

        Parameters
        ----------
        i : int
            Row position to insert the markers into the sample alignment.
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        markers : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if not(isinstance(ids, list) and
               sum((isinstance(j, str) for j in ids))):
            raise TypeError('ids must be a list of str.')
        if not(isinstance(descriptions, list) and
               sum((isinstance(j, str) for j in descriptions))):
            raise TypeError('descriptions must be a list of str.')
        if not(isinstance(markers, list) and
               sum((isinstance(j, str) for j in markers))):
            raise TypeError('markers must be a list of str.')
        aln.markers.insert_rows(i, ids, descriptions, markers)
        if copy:
            return aln

    def append_markers_from_lists(self, ids, descriptions, markers, copy=False):
        """Inserts new sequences after the last row of the alignment matrix
        inplace. This increases the total number of markers.

        Parameters
        ----------
        ids : list of str
            List of new sample names/IDs to be appended to the
            existing list of sample names.
        descriptions : list of str
            List of new sample descriptions to be appended to the
            existing list of descriptions. Length and order must correspond to
            the given list of ids.
        markers : list of str
            List of new sample sequences to be appended to the
            existing list of sequences. Length and order must correspond to
            the given list of ids.
        copy : bool, optional
            Returns a new copy instead of appending the markers inplace.
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if not(isinstance(ids, list) and
               sum((isinstance(j, str) for j in ids))):
            raise TypeError('ids must be a list of str.')
        if not(isinstance(descriptions, list) and
               sum((isinstance(j, str) for j in descriptions))):
            raise TypeError('descriptions must be a list of str.')
        if not(isinstance(markers, list) and
               sum((isinstance(j, str) for j in markers))):
            raise TypeError('markers must be a list of str.')
        aln.markers.append_rows(ids, descriptions, markers)
        if copy:
            return aln


    # Setter methods
    # ==========================================================================
    def replace_samples(self, i, sequences, copy=False):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the samples to be replaced.
            If a list, length must match the length of `sequences`.
        sequences : str
            str or list of new sequences. If `i` is an int or str,
            must be a str. If `i` is a list, must be a list and
            length must match the length of `i`.
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if isinstance(i, int) and isinstance(sequences, str):
            aln.samples.set_sequences([i], [sequences])
        elif isinstance(i, str) and isinstance(sequences, str):
            ids = aln.samples.row_names_to_ids([i])
            aln.samples.set_sequences([ids], [sequences])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.samples.set_sequences(i, sequences)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            ids = aln.samples.row_names_to_ids(i)
            aln.samples.set_sequences(ids, sequences)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if copy:
            return aln

    def replace_markers(self, i, sequences, copy=False):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the markers to be replaced.
            If a list, length must match the length of `sequences`.
        sequences : str
            str or list of new sequences. If `i` is an int or str,
            must be a str. If `i` is a list, must be a list and
            length must match the length of `i`.
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Calls specific set_sequence setter depending on the
        # type if i
        if isinstance(i, int) and isinstance(sequences, str):
            aln.markers.set_sequences([i], [sequences])
        elif isinstance(i, str) and isinstance(sequences, str):
            ids = aln.markers.row_names_to_ids([i])
            aln.markers.set_sequences([ids], [sequences])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.markers.set_sequences(i, sequences)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            ids = aln.markers.row_names_to_ids(i)
            aln.markers.set_sequences(ids, sequences)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if copy:
            return aln

    def reorder_samples(self, ids, copy=False):
        """Reorders samples based on a list of names or positions.

        Parameters
        ----------
        ids : list of int or list of str
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Check type of i, and convert if necessary
        if isinstance(ids, list) and sum((isinstance(j, int) for j in ids)):
            pass
        elif isinstance(ids, list) and sum((isinstance(j, str) for j in ids)):
            ids = aln.samples.row_names_to_ids(ids)
        else:
            raise TypeError('ids must be a list of int or list of str.')
        aln.samples.reorder_rows(ids)
        if copy:
            return aln

    def reorder_markers(self, ids, copy=False):
        """Reorders markers based on a list of names or positions.

        Parameters
        ----------
        ids : list of int or list of str
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Check type of i, and convert if necessary
        if isinstance(ids, list) and sum((isinstance(j, int) for j in ids)):
            pass
        elif isinstance(ids, list) and sum((isinstance(j, str) for j in ids)):
            ids = aln.markers.row_names_to_ids(ids)
        else:
            raise TypeError('ids must be a list of int or list of str.')
        aln.markers.reorder_rows(ids)
        if copy:
            return aln

    def reset_coordinates(self, start=0, stop=None, state=1):
        """Resets the coordinates associated to the alignment columns.  

        Parameters
        ----------
        start : int, optional
            Starting coordinate value (default is 0)
        stop : int, optional
            End coordinate value. This value is not part of linear space
            represented by the alignment length. (default is None, ending value
            becomes the sum of start and the number of sites)
        state : int, optional
            State annotation. (default is 1)

        """
        start = start if start is not None else 0
        stop = stop if stop is not None else start + self.samples.nsites
        state = state if state is not None else 1
        self._linspace: BlockSpace = BlockSpace(start, stop, state)


    # Deleter methods
    # ==========================================================================

    # General deleters
    # ------------------------------
    def remove_sites(self, i, copy=False):
        """Removes sites based on the given list of column numbers.

        This is the functional opposite of the `retain_sites` method.

        Parameters
        ----------
        i : int, or list of int
            int or list specifying the sites to be removed.
        copy : bool, optional
            Returns a new copy instead of removing sites inplace.
            (default is False, operation is done inplace)

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Check type of i, and convert if necessary
        if isinstance(i, int):
            i = [i]
        # Perform removal inplace
        aln.samples.remove_sites(i)
        if aln.markers:
            aln.markers.remove_sites(i)
            assert aln.samples.nsites == aln.markers.nsites, \
                "Sample and marker nsites are not equal."
        aln._linspace.remove(i)
        if copy:
            return aln

    def retain_sites(self, i, copy=False):
        """Keeps sites based on the given list of column numbers.

        This is the functional opposite of the `remove_sites` method.

        Parameters
        ----------
        i : int or list of int
            int or list specifying the sites to be retained.
        copy : bool, optional
            Returns a new copy instead of performing the operation inplace.
            (default is False, operation is done inplace)

        Returns
        -------
        Alignment or None
            If copy is True, returns a new alignment, otherwise no
            value is returned (None).

        """
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        # Check type of i, and convert if necessary
        if isinstance(i, int):
            i = [i]
        # Perform removal inplace
        aln.samples.retain_sites(i)
        if aln.markers:
            aln.markers.retain_sites(i)
            assert aln.samples.nsites == aln.markers.nsites, \
                "Sample and marker nsites are not equal."
        aln._linspace.retain(i)
        if copy:
            return aln

    # Sample deleters
    # ------------------------------
    def remove_samples(self, i, match_prefix=False, match_suffix=False,
                       copy=False):
        """Removes sample sequences based on the given index.

        This is the functional opposite of the `retain_samples` method.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the samples to be removed.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        if isinstance(i, int):
            aln.samples.remove_rows([i])
        elif isinstance(i, str):
            if match_prefix:
                aln.samples.remove_rows_by_prefix([i])
            elif match_suffix:
                aln.samples.remove_rows_by_suffix([i])
            else:
                aln.samples.remove_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.samples.remove_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                aln.samples.remove_rows_by_prefix(i)
            elif match_suffix:
                aln.samples.remove_rows_by_suffix(i)
            else:
                aln.samples.remove_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if copy:
            return aln

    def retain_samples(self, i, match_prefix=False, match_suffix=False, copy=False):
        """Keeps sample sequences based on the given index.

        This is the functional opposite of the `remove_samples` method.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the samples to be retained.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
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
        aln = self.__class__(
            self.name, self.samples.copy(), self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        if isinstance(i, int):
            aln.samples.retain_rows([i])
        if isinstance(i, str):
            if match_prefix:
                aln.samples.retain_rows_by_prefix([i])
            elif match_suffix:
                aln.samples.retain_rows_by_suffix([i])
            else:
                aln.samples.retain_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.samples.retain_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                aln.samples.retain_rows_by_prefix(i)
            elif match_suffix:
                aln.samples.retain_rows_by_suffix(i)
            else:
                aln.samples.retain_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if copy:
            return aln

    # Marker deleters
    # ------------------------------
    def remove_markers(self, i, match_prefix=False, match_suffix=False,
                       copy=False):
        """Removes sample sequences based on the given index.

        This is the functional opposite of the `retain_markers` method.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the markers to be removed.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
            the list of sample names. This parameter is considered
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        if isinstance(i, int):
            aln.markers.remove_rows([i])
        elif isinstance(i, str):
            if match_prefix:
                aln.markers.remove_rows_by_prefix([i])
            elif match_suffix:
                aln.markers.remove_rows_by_suffix([i])
            else:
                aln.markers.remove_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.markers.remove_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                aln.markers.remove_rows_by_prefix(i)
            elif match_suffix:
                aln.markers.remove_rows_by_suffix(i)
            else:
                aln.markers.remove_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
        if copy:
            return aln

    def retain_markers(self, i, match_prefix=False, match_suffix=False, copy=False):
        """Keeps sample sequences based on the given index.

        This is the functional opposite of the `remove_markers` method.

        Parameters
        ----------
        i : int, str, list of int, or list of str
            int, str or list specifying the markers to be retained.
        match_prefix : bool, optional
            Whether to interpret `i` as a prefix to match against
            the list of sample names. (default is False)
        match_suffix : bool, optional
            Whether to interpret `i` as a suffix to match against
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
        aln = self.__class__(
            self.name,
            self.samples.copy(),
            self.markers.copy(),
            _linspace=self._linspace.copy()) if copy else \
            self
        if isinstance(i, int):
            aln.markers.retain_rows([i])
        if isinstance(i, str):
            if match_prefix:
                aln.markers.retain_rows_by_prefix([i])
            elif match_suffix:
                aln.markers.retain_rows_by_suffix([i])
            else:
                aln.markers.retain_rows_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            aln.markers.retain_rows(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                aln.markers.retain_rows_by_prefix(i)
            elif match_suffix:
                aln.markers.retain_rows_by_suffix(i)
            else:
                aln.markers.retain_rows_by_name(i)
        else:
            raise TypeError('i must be an int, str, list of int, or list of str.')
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
            Stopping column position. If None, the iterator will continue
            until the end of the alignment. (default is None)
        size : int, optional
            Size of chunks to yield. For single characters, `size` = 1.
            For codons, `size` = 3. (default is 1)

        Raises
        ------
        ValueError
            If the alignment cannot be cleanly cut up into the specified
            chunk size (`nsites` not divisible be `size`),
            a ValueError is raised.

        Yields
        ------
        list of str
            List of sequences representing a site or chunk of the alignment.

        """
        if stop is None:
            stop = self.nsites
        if (stop - start) % size != 0:
            raise ValueError('Alignment cannot be completely divided into '
                             'chucks of size {}'.format(size))
        for i in range(start, stop, size):
            samples = [s[i:i+size] for s in self.samples.sequences]
            if not (self.markers is None or self.markers.nrows == 0):
                markers = [s[i:i+size] for s in self.markers.sequences]
                yield samples + markers
            else:
                yield samples

    def iter_sample_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the sample alignment. Excludes markers.

        Parameters
        ----------
        start : int, optional
            Starting position. (default is 0)
        stop : [type], optional
            Stopping column position. If None, the iterator will continue
            until the end of the alignment. (default is None)
        size : int, optional
            Size of chunks to yield. For single characters, `size` = 1.
            For codons, `size` = 3. (default is 1)

        Raises
        ------
        ValueError
            If the alignment cannot be cleanly cut up into the specified
            chunk size (`nsites` not divisible be `size`),
            a ValueError is raised.

        Yields
        ------
        list of str
            List of sequences representing a site or chunk of the alignment.

        """
        if stop is None:
            stop = self.nsites
        if (stop - start) % size != 0:
            raise ValueError('Alignment cannot be completely divided into '
                             'chucks of size {}'.format(size))
        for i in range(start, stop, size):
            yield [s[i:i+size] for s in self.samples.sequences]

    def iter_marker_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the marker alignment. Excludes samples.

        Parameters
        ----------
        start : int, optional
            Starting position. (default is 0)
        stop : [type], optional
            Stopping column position. If None, the iterator will continue
            until the end of the alignment. (default is None)
        size : int, optional
            Size of chunks to yield. For single characters, `size` = 1.
            For codons, `size` = 3. (default is 1)

        Raises
        ------
        ValueError
            If the alignment cannot be cleanly cut up into the specified
            chunk size (`nsites` not divisible be `size`),
            a ValueError is raised.

        Yields
        ------
        list of str
            List of sequences representing a site or chunk of the alignment.

        """
        if self.markers is None or self.markers.nrows == 0:
            return
        if stop is None:
            stop = self.nsites
        if (stop - start) % size != 0:
            raise ValueError('Alignment cannot be completely divided into '
                             'chucks of size {}'.format(size))
        for i in range(start, stop, size):
            yield [s[i:i+size] for s in self.markers.sequences]


    # Format converters
    # ==========================================================================
    @classmethod
    def from_fasta(cls, path, name=None, marker_kw=None, comment_parser=None):
        """Create an Alignment object from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file.
        name : str
            Name of the new alignment.
        marker_kw : str, optional
            A sample is considered a marker if this keyword is found
            in the identifier.

        Returns
        -------
        Alignment
            Creates a new Alignment object based on the identifiers,
            descriptions, and sequences in the FASTA file.

        """
        if marker_kw is None:
            marker_kw = ''
        # Create alignments
        samples, markers, comment_list = \
            fasta_file_to_basealignments(path, marker_kw)
        comment_parser = \
            comment_parser if comment_parser is not None else parse_comment_list
        kwargs = comment_parser(comment_list)
        name = name if name else kwargs['name']
        return cls(name, samples, markers, **kwargs)

    def to_fasta(self, path, include_markers=True):
        """Saves the alignment as a FASTA-formatted text file.

        Parameters
        ----------
        path : str
            Path to save the alignment to.

        """
        with open(path, 'w') as writer:
            print(';name\t' + str(self.name), file=writer)
            print(';coords\t{' + self._linspace.to_simple_block_str() + '}',
                  file=writer)
            print(self.samples, file=writer)
            if include_markers:
                print(self.markers, file=writer)


    # Special methods
    # ==========================================================================
    def __getitem__(self, key):
        if key in self.samples.ids():
            i = self.samples.row_names_to_ids([key])[0]
            return self.samples.get_row(i)
        elif key in self.markers.ids():
            i = self.markers.row_names_to_ids([key])[0]
            return self.markers.get_row(i)
        raise KeyError('Key did not match any sample or marker ID')

    def __delitem__(self, key):
        if key in self.samples.ids():
            i = self.samples.row_names_to_ids([key])
            return self.samples.remove_rows(i)
        elif key in self.markers.ids():
            i = self.markers.row_names_to_ids([key])
            return self.markers.remove_rows(i)
        raise KeyError('Key did not match any sample or marker ID')

    def __iter__(self):
        yield from self.iter_sites()

    def __repr__(self):
        return '{}(nsamples={}, nsites={}, nmarkers={})'.format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            self.nmarkers
        )

    def __str__(self):
        if self.markers:
            return '\n'.join([
                ';name\t' + str(self.name),
                ';coords\t{' + self._linspace.to_simple_block_str() + '}',
                str(self.samples),
                str(self.markers),
            ])
        return '\n'.join([
            ';name\t' + str(self.name),
            ';coords\t{' + self._linspace.to_simple_block_str() + '}',
            str(self.samples),
        ])

    def __len__(self):
        raise NotImplementedError(
            'len() is not implemented for Alignment.\n'
            'Use .nsites to get the number of columns, '
            '.nsamples to get the number of samples, '
            '.nmarkers to get the number of markers, or '
            '.nrows to get all the number of alignment rows.')
