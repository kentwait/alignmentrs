from collections.abc import Iterable
from collections import namedtuple
from copy import deepcopy
import itertools
import numpy as np
from alignmentrs import Sequence, Marker
from blockrs import Block
from blockrs.block import remove_sites

vector_ord = np.vectorize(ord)
vector_chr = np.vectorize(chr)

class AlignmentMatrix(object):
    """AlignmentMatrix depicts an aligned set of sequences as
    a 2d array of uint32 values.
    """
    def __init__(self, nsamples, nsites, to_uint_fn=None, from_uint_fn=None):
        """Creates a new alignment matrix.

        Parameters
        ----------
        nsamples : int
            Number of samples (rows)
        nsites : int
            Number of sites (columns)

        """
        self.custom_to_uint_fn = None if to_uint_fn is None else \
                                 np.vectorize(to_uint_fn)
        self.custom_from_uint_fn = None if from_uint_fn is None else \
                                   np.vectorize(from_uint_fn)
        self.matrix = np.empty((nsamples, nsites), dtype=np.uint32)

    def to_uint(self, x):
        if self.custom_to_uint_fn is None:
            return vector_ord(x)
        return self.custom_to_uint_fn(x)

    def from_uint(self, x):
        if self.custom_from_uint_fn is None:
            return vector_chr(x)
        return self.custom_from_uint_fn(x)

    @classmethod
    def empty(cls, nsamples, nsites):
        """Creates an empty alignment matrix with a given number of rows
        (nsamples) and columns (nsites).

        Parameters
        ----------
        nsamples : int
            Number of samples (rows)
        nsites : int
            Number of sites (columns)

        """
        return cls(nsamples, nsites)

    @property
    def nsamples(self):
        """Returns the number of samples (rows) in the alignment matrix.
        """
        return self.matrix.shape[0]

    @property
    def nsites(self):
        """Returns the number of sites (columns) in the alignment matrix.
        """
        return self.matrix.shape[1]

    @property
    def shape(self):
        """Returns the shape of the alignment matrix.
        
        Returns
        -------
        tuple
            Tuple of number of rows and columns of the alignment

        """
        return self.matrix.shape

    @classmethod
    def subset(cls, aln_matrix, rows=None, cols=None,
               row_step=1, col_step=1):
        """Returns a subset of the alignment matrix by both samples and sites.

        Parameters
        ----------
        m : AlignmentMatrix
        rows : list
        cols : list
        row_step : int
        col_step : int

        Returns
        -------
        AlignmentMatrix

        """
        if rows is None:
            rows = range(0, aln_matrix.nsamples, row_step)
        else:
            if isinstance(rows, int):
                rows = [rows]
            if row_step != 1:
                raise ValueError('row_step value is considered only if rows ' \
                                 'is None')
        if cols is None:
            cols = range(0, aln_matrix.nsites, col_step)
        else:
            if isinstance(cols, int):
                cols = [cols]
            if col_step != 1:
                raise ValueError('col_step value is considered only if cols ' \
                                 'is None')
        new_matrix = np.copy(aln_matrix.matrix[rows][:, cols])
        # Create a new alignment matrix
        new_aln_matrix = cls(*new_matrix.shape)
        # Overwrite contents
        new_aln_matrix.custom_to_uint_fn = deepcopy(aln_matrix.to_uint_fn)
        new_aln_matrix.custom_from_uint_fn = deepcopy(aln_matrix.from_uint_fn)
        new_aln_matrix.matrix = new_matrix
        return new_aln_matrix

    @classmethod
    def from_uint_matrix(cls, matrix, to_uint_fn=None, from_uint_fn=None):
        """Creates a new AlignmentMatrix from a numpy array.

        Parameters
        ----------
        matrix : np.array
            Numpy array containg the data
        to_uint_fn : function, optional
            Function to convert data into uint32
        from_uint_fn : [type], optional
            Function to convert uint32 encoding back to the original data

        Returns
        -------
        AlignmentMatrix

        """
        # Create a new alignment matrix
        new_matrix = cls(*matrix.shape, to_uint_fn=to_uint_fn,
                         from_uint_fn=from_uint_fn)
        # Overwrite and copy uint_matrix
        new_matrix.matrix = np.copy(matrix)
        return new_matrix

    def replace_sample(self, i, sequence):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        sequence : str or iterable
        i : int

        """
        # Check if nsites is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsites))
        if isinstance(sequence, str):
            self.matrix[i, :] = self.to_uint(list(sequence))
        elif isinstance(sequence, Iterable):
            self.matrix[i, :] = self.to_uint(sequence)
        else:
            raise TypeError('sequence must be a string or an iterable')

    def insert_sample(self, i, sequence):
        """Inserts a new sequence in the alignment matrix at the specified
        row position. This increases the total number of rows.

        Parameters
        ----------
        sequence : str or iterable
        i : int

        """
        # Check if nsites is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsites))
        if isinstance(sequence, str):
            new_array = np.array([self.to_uint(list(sequence))],
                                 dtype=np.uint32)
        elif isinstance(sequence, Iterable):
            new_array = np.array([self.to_uint(sequence)], dtype=np.uint32)
        else:
            raise TypeError('sequence must be a string or an iterable')
        self.matrix = np.insert(self.matrix, i, new_array, axis=0)

    def append_sample(self, sequence):
        """Inserts a new sequence after the last row of the alignment matrix.
        This increases the total number of rows by 1.

        Parameters
        ----------
        sequence : str or iterable

        """
        # Check if nsites is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsites))
        if isinstance(sequence, str):
            new_array = np.array([self.to_uint(list(sequence))],
                                 dtype=np.uint32)
        elif isinstance(sequence, Iterable):
            new_array = np.array([self.to_uint(sequence)], dtype=np.uint32)
        else:
            raise TypeError('sequence must be a string or an iterable')
        self.matrix = np.append(self.matrix, new_array, axis=0)

    def remove_samples(self, i):
        """Removes sample sequences based on the given index.
        If index is a number, only one sequence is removed.
        If the index is a list of numbers, the sequence found at each row
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        self.matrix = np.delete(self.matrix, i, axis=0)

    def get_samples(self, i):
        """Returns a new alignment matrix containing only the samples specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        AlignmentMatrix

        """
        if isinstance(i, int):
            i = [i]
        return self.__class__.subset(self, rows=i)

    def get_samples_as_str(self, i):
        """Returns a list of sequence strings containing only the samples
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        if isinstance(i, int):
            i = [i]
        translated_g = (self.from_uint(row) for row in self.matrix[i, :])
        return [''.join(trans) for trans in translated_g]

    def get_sites(self, i):
        """Returns a new alignment matrix containing only the sites specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        AlignmentMatrix

        """
        if isinstance(i, int):
            i = [i]
        return self.__class__.subset(self, cols=i)

    def get_sites_as_str(self, i):
        """Returns a list of sequence strings containing only the sites
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        if isinstance(i, int):
            i = [i]
        translated_g = (self.from_uint(row) for row in self.matrix[:, i])
        return [''.join(trans) for trans in translated_g]

    def replace_site(self, i, sequence):
        """Replaces the sequence for a given column in the alignment matrix.

        Parameters
        ----------
        sequence : str or iterable
        i : int

        """
        # Check if nsamples is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsamples))
        if isinstance(sequence, str):
            self.matrix[:, i] = self.to_uint(list(sequence))
        elif isinstance(sequence, Iterable):
            self.matrix[:, i] = self.to_uint(sequence)
        else:
            raise TypeError('sequence must be a string or an iterable')

    def insert_site(self, i, sequence):
        """Inserts a new sequence in the alignment matrix at the specified
        site position. This increases the total number of columns.

        Parameters
        ----------
        sequence : str or iterable
        i : int

        """
        # Check if nsamples is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsamples))
        if isinstance(sequence, str):
            new_array = np.array([self.to_uint(list(sequence))],
                                 dtype=np.uint32)
        elif isinstance(sequence, Iterable):
            new_array = np.array([self.to_uint(sequence)], dtype=np.uint32)
        else:
            raise TypeError('sequence must be a string or an iterable')
        self.matrix = np.insert(self.matrix, i, new_array, axis=1)

    def append_site(self, sequence):
        """Inserts a new sequence at after the last column of the
        alignment matrix. This increases the total number of columns by 1.

        Parameters
        ----------
        sequence : str or iterable

        """
        # Check if nsamples is equal
        if len(sequence) != self.nsites:
            raise ValueError(
                'length of sequence not equal to {}'.format(self.nsamples))
        if isinstance(sequence, str):
            new_array = np.array([self.to_uint(list(sequence))],
                                 dtype=np.uint32).T
        elif isinstance(sequence, Iterable):
            new_array = np.array([self.to_uint(sequence)], dtype=np.uint32).T
        else:
            raise TypeError('sequence must be a string or an iterable')
        self.matrix = np.append(self.matrix, new_array, axis=1)

    def remove_sites(self, i):
        """Removes sites based on the given index.
        If index is a number, only one site is removed.
        If the index is a list of numbers, the sequence found at each column
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        self.matrix = np.delete(self.matrix, i, axis=1)

    def to_string_list(self):
        return [''.join(self.from_uint(row)) for row in self.matrix]

    def __len__(self):
        return len(self.matrix)

    def __repr__(self):
        return "{}(nsamples={}, nsites={})\n{}".format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            repr(self.from_uint(self.matrix))
        )

    def __str__(self):
        return '\n'.join(self.get_samples_as_str(list(range(len(self)))))


class BaseAlignment(AlignmentMatrix):
    """BaseAlignment represents a mulitple sequence alignment
    by storing the IDs and descriptions in individual lists, and
    representing the aligned sequences as a 2d array of uint32 values.
    """
    def __init__(self, sequence_list, to_uint_fn=None, from_uint_fn=None):
        """Creates a new alignment from a list of Sequence objects

        Parameters
        ----------
        sequence_list : list of Sequence
        to_uint_fn : function
        from_uint_fn : function

        """
        nsamples = len(sequence_list)
        nsites = len(sequence_list[0].sequence) if nsamples > 0 else 0
        # Call parent constructor
        super().__init__(nsamples, nsites, to_uint_fn=to_uint_fn,
                         from_uint_fn=from_uint_fn)
        # Add new attributes
        self.ids = [s.id for s in sequence_list]
        self.descriptions = [s.description for s in sequence_list]
        self._metadata = ('ids', 'descriptions')
        # Overwrite empty matrix with data
        self.matrix = np.array([self.to_uint(list(s.sequence))
                                for s in sequence_list], dtype=np.uint32)

    @classmethod
    def subset(cls, aln, rows=None, cols=None, row_step=1, col_step=1):
        """Returns a subset of the alignment matrix by both samples and sites.

        Parameters
        ----------
        m : BaseAlignment
        rows : list
        cols : list
        row_step : int
        col_step : int

        Returns
        -------
        BaseAlignment

        """
        if rows is None:
            rows = range(0, aln.nsamples, row_step)
        else:
            if isinstance(rows, int):
                rows = [rows]
            if row_step != 1:
                raise ValueError('row_step value is considered only if rows' \
                                 'is None')
        if cols is None:
            cols = range(0, aln.nsites, col_step)
        else:
            if isinstance(cols, int):
                cols = [cols]
            if col_step != 1:
                raise ValueError('col_step value is considered only if cols ' \
                                 'is None')
        # Create a new BaseAlignment
        new_aln = cls.__new__(cls)
        # Add attributes
        new_aln.custom_from_uint_fn = deepcopy(aln.custom_from_uint_fn)
        new_aln.custom_to_uint_fn = deepcopy(aln.custom_to_uint_fn)
        new_aln.matrix = np.copy(aln.matrix[rows][:, cols])
        # Add metadata
        new_aln._metadata = aln._metadata  # copy tuple
        for name in new_aln._metadata:
            new_value = [v for i, v in enumerate(aln.__getattribute__(name))
                         if i in rows]
            new_aln.__setattr__(name, new_value)
        return new_aln

    def get_samples(self, i):
        """Returns a new alignment matrix containing only the samples specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        BaseAlignment

        """
        if isinstance(i, int):
            return self.__class__.subset(self, rows=i)
        # handles single string input
        elif isinstance(i,str):
            i = self.ids.index(i) # gets id index from id
            return self.__class__.subset(self, rows=i)  #retrieve info with id
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            elif len(i) == sum((isinstance(j, str) for j in i)):
                # get int position in self.ids and then use subset
                i = [self.ids.index(idx) for idx in i]
            else:
                raise ValueError('i must be a list of int or str')
        else:
            raise ValueError('i must be an int or str or a list of int or str')
        return self.__class__.subset(self, rows=i)

    def get_samples_as_str(self, i):
        """Returns a list of sequences as strings containing only the
        specified samples.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        sample_array = self.get_samples(i)
        translated_g = (self.from_uint(row) for row in sample_array.matrix)
        return [''.join(trans) for trans in translated_g]

    def remove_samples(self, i):
        """Removes sample sequences based on the given index.
        If index is a number, only one sequence is removed.
        If the index is a list of numbers, the sequence found at each row
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        if isinstance(i, int):
            return self.__class__.subset(self, rows=i)
        # same as above
        elif isinstance(i,str):
            i = self.ids.index(i) # gets id index from id
            return self.__class__.subset(self, rows=i)  # retrieve info with id
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            elif len(i) == sum((isinstance(j, str) for j in i)):
                # get int position in self.ids and then use subset
                i = [self.ids.index(idx) for idx in i]
            else:
                raise ValueError('i must be a list of int or str')
        else:
            raise ValueError('i must be an int or str or a list of int or str')
        super().remove_samples(i)

    def get_sites(self, i):
        """Returns a new alignment matrix containing only the sites specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        AlignmentMatrix

        """
        if isinstance(i, int):
            return self.__class__.subset(self, rows=i)
        # same as above
        elif isinstance(i,str):
            i = self.ids.index(i) # gets id index from id
            return self.__class__.subset(self, rows=i)  # retrieve info with id
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            else:
                raise ValueError('i must be a list of int')
        else:
            raise ValueError('i must be an int or str or a list of int or str')
        return self.__class__.subset(self, cols=i)

    def get_sites_as_str(self, i):
        """Returns a list of sequences as strings containing only the
        specified sites.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        sample_array = self.get_sites(i)
        translated_g = (self.from_uint(row) for row in sample_array.matrix)
        return [''.join(trans) for trans in translated_g]

    def remove_sites(self, i):
        """Removes sites based on the given index.
        If index is a number, only one site is removed.
        If the index is a list of numbers, the sequence found at each column
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        super().remove_sites(i)

    @classmethod
    def from_uint_matrix(cls, matrix, ids, descriptions,
                         to_uint_fn=None, from_uint_fn=None):
        # Create an empty BaseAlignment
        new_aln = cls.__new__(cls)
        # Assign uint conversion functions
        new_aln.custom_from_uint_fn = from_uint_fn
        new_aln.custom_to_uint_fn = to_uint_fn
        # Assign values
        new_aln.ids = deepcopy(ids)
        new_aln.descriptions = deepcopy(descriptions)
        new_aln._metadata = ('ids', 'descriptions')
        new_aln.matrix = np.copy(matrix)
        return new_aln

    def __str__(self):
        id_lines = ['>{} {}'.format(sid, desc) if desc else '>{}'.format(sid)
                    for sid, desc in zip(self.ids, self.descriptions)]
        seq_lines = self.get_samples_as_str(list(range(len(self))))
        return '\n'.join(itertools.chain(*zip(id_lines, seq_lines)))


class SampleAlignment(BaseAlignment):
    """SampleAlignment represents a mulitple sequence alignment of
    biological sequences.
    """
    def __init__(self, sequence_list, to_uint_fn=None, from_uint_fn=None,
                 to_block_fn=None, from_block_fn=None):
        super().__init__(sequence_list, to_uint_fn=to_uint_fn,
                         from_uint_fn=from_uint_fn)
        # Set conversion functions
        self.custom_to_block_fn = to_block_fn
        self.custom_from_block_fn = from_block_fn
        # Generate sample block lists from sample descriptions
        self.block_lists = [self.to_block(desc)
                            for desc in self.descriptions]
        # Include block_lists in _metadata
        self._metadata = ('ids', 'descriptions', 'block_lists')

    def to_block(self, string):
        if self.custom_to_block_fn is None:
            tuple_list = (tuple(map(int, paired.split(':')))
                          for paired in string.split('_')[-1].split(';'))
            return [Block(tpl[0], tpl[1]) for tpl in tuple_list]
        return self.custom_to_block_fn(string)

    def from_block(self, block_list):
        if self.custom_from_block_fn is None:
            return '{}_{}'.format(
                len(block_list), ';'.join([str(b) for b in block_list]),
            )
        return self.custom_from_block_fn(block_list)

    @classmethod
    def subset(cls, aln, rows=None, cols=None, row_step=1, col_step=1):
        if rows is None:
            rows = range(0, aln.nsamples, row_step)
        else:
            if isinstance(rows, int):
                rows = [rows]
            if row_step != 1:
                raise ValueError('row_step value is considered only if rows' \
                                 'is None')
        if cols is None:
            cols = range(0, aln.nsites, col_step)
        else:
            if isinstance(cols, int):
                cols = [cols]
            if col_step != 1:
                raise ValueError('col_step value is considered only if cols ' \
                                 'is None')
        def f(seq, blist, drop_pos_lst):
            _, new_blist = remove_sites(seq, blist, drop_pos_lst)
            return new_blist

        # Creates a new SampleAlignment
        new_aln = cls.__new__(cls)
        # Copies custom functions
        new_aln.custom_from_uint_fn = deepcopy(aln.custom_from_uint_fn)
        new_aln.custom_to_uint_fn = deepcopy(aln.custom_to_uint_fn)
        new_aln.custom_to_block_fn = deepcopy(aln.custom_to_block_fn)
        new_aln.custom_from_block_fn = deepcopy(aln.custom_from_block_fn)
        # Copy the subset of the matrix
        new_aln.matrix = np.copy(aln.matrix[rows][:, cols])
        # Copies metadata
        new_aln.ids = [v for i, v in enumerate(aln.ids) if i in rows]
        new_aln._metadata = aln._metadata  # copy tuple
        # Update block_lists
        drop_pos_lst = [i for i in range(aln.nsites) if i not in cols]
        seq_list = (row for row in aln.to_string_list())
        new_aln.block_lists = [
            f(seq, blist, drop_pos_lst)
            for seq, blist in zip(seq_list, aln.block_lists)
        ]
        # Update descriptions
        new_aln.descriptions = [new_aln.from_block(blist)
                                for blist in new_aln.block_lists]
        return new_aln

    @classmethod
    def from_uint_matrix(cls, matrix, ids, descriptions, block_lists,
                         to_uint_fn=None, from_uint_fn=None,
                         to_block_fn=None, from_block_fn=None):
        # Create an empty SampleAlignment
        new_aln = cls.__new__(cls)
        # Assign uint conversion functions
        new_aln.custom_to_uint_fn = None if to_uint_fn is None else \
                                    np.vectorize(to_uint_fn)
        new_aln.custom_from_uint_fn = None if from_uint_fn is None else \
                                    np.vectorize(from_uint_fn)
        new_aln.custom_to_block_fn = to_block_fn
        new_aln.custom_from_block_fn = from_block_fn
        # Assign values
        new_aln.ids = deepcopy(ids)
        new_aln.descriptions = deepcopy(descriptions)
        new_aln.block_lists = [[Block(b.start, b.stop) for b in blist]
                               for blist in block_lists]  # deepcopy blocks
        new_aln._metadata = ('ids', 'descriptions', 'block_lists')
        # Copy the matrix
        new_aln.matrix = np.copy(matrix)
        return new_aln

    # TODO: Add a from_fasta method

    # TODO: Add a method that adds markers and makes this into an Alignment object


class MarkerAlignment(BaseAlignment):
    """MarkerAlignment represents one or more alignment markers of
    a multiple sequence alignment.
    """
    @property
    def nmarkers(self):
        """Returns the number of markers in the alignment.
        """
        return self.nsamples

    def get_markers(self, i):
        """Returns a new alignment matrix containing only the markers specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        MarkerAlignment

        """
        return self.get_samples(i)

    def get_markers_as_str(self, i):
        """Returns a list of sequences as strings containing only the
        specified markers.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        return self.get_samples_as_str(i)

    def remove_markers(self, i):
        """Removes sites based on the given index.
        If index is a number, only one site is removed.
        If the index is a list of numbers, the sequence found at each column
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        return self.remove_samples(i)


class BinaryMarkerAlignment(MarkerAlignment):
    def __init__(self, marker_list):
        """Creates a BinaryMarkerAlignment from a list of of Marker

        Parameters
        ----------
        marker_list : list of Marker

        """
        to_uint_fn = int  # assumes x is '0' or '1'
        from_uint_fn = str   # assumes x is int 0 or 1
        super().__init__(marker_list, to_uint_fn=to_uint_fn, 
                         from_uint_fn=from_uint_fn)


class Alignment(object):
    """Alignment is a complete representation of a multiple sequence alignment
    of biological sequences an their annotations such as alignment markers and
    alignment block data.
    """
    def __init__(self, sequence_list, marker_list, name=None,
                 sample_to_uint_fn=None, uint_to_sample_fn=None,
                 marker_to_uint_fn=None, uint_to_marker_fn=None):
        """Creates a new Alignment object from a list of Sequence and Marker objects.

        Parameters
        ----------
        sequence_list : list of Sequence
        marker_list : list of Marker

        """
        self.name = name
        self._sample_aln = SampleAlignment(sequence_list,
                                           to_uint_fn=sample_to_uint_fn,
                                           from_uint_fn=uint_to_sample_fn)
        self._marker_aln = MarkerAlignment(marker_list,
                                           to_uint_fn=marker_to_uint_fn,
                                           from_uint_fn=uint_to_marker_fn)
        # assert self._sample_aln.nsites == self._marker_aln.nsites

    @property
    def nsites(self):
        """Returns the number of sites in the alignment.
        """
        return self._sample_aln.nsites

    @property
    def nsamples(self):
        """Returns the number of samples in the alignment.
        """        
        return self._sample_aln.nsamples

    @property
    def nmarkers(self):
        """Returns the number of markers in the alignment
        """
        return self._marker_aln.nmarkers

    @property
    def samples(self):
        """Returns the sample alignment
        """
        return self._sample_aln

    @property
    def sample_matrix(self):
        """Returns the sample alignent matrix
        """
        return self._sample_aln.matrix

    @property
    def markers(self):
        """Returns the marker alignment
        """
        return self._marker_aln

    @property
    def marker_matrix(self):
        """Returns the marker alignment matrix
        """
        return self._marker_aln.matrix

    @classmethod
    def subset(cls, aln, sample_ids=None, marker_ids=None, sites=None,
               sample_id_step=1, marker_id_step=1, site_step=1):
        """Returns a subset of the alignment by samples, markers and sites.

        Parameters
        ----------
        aln : Alignment
        sample_ids : list
        marker_ids : list
        sites : list
        sample_id_step : int
        marker_id_step : int
        site_step : int

        Returns
        -------
        Alignment

        """
        if sample_ids is None:
            sample_ids = range(0, aln.nsamples, sample_id_step)
        else:
            if sample_id_step != 1:
                raise ValueError('sample_id_step value is considered only ' \
                                 'if sample_ids is None')
        if marker_ids is None:
            marker_ids = range(0, aln.nmarkers, marker_id_step)
        else:
            if marker_id_step != 1:
                raise ValueError('marker_id_step value is considered only ' \
                                 'if marker_ids is None')
        if sites is None:
            sites = range(0, aln.nsites, site_step)
        else:
            if site_step != 1:
                raise ValueError('site_step value is considered only ' \
                                 'if sites is None')
        new_aln = cls.__new__(cls)
        new_aln.name = aln.name
        new_aln._sample_aln = aln._sample_aln.__class__.subset(
            aln._sample_aln,
            rows=sample_ids, cols=sites,
            row_step=sample_id_step, col_step=site_step
        )
        new_aln._marker_aln = aln._marker_aln.__class__.subset(
            aln._marker_aln,
            rows=marker_ids, cols=sites,
            row_step=marker_id_step, col_step=site_step
        )
        return new_aln

    def replace_sample(self, sequence_str, i):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        sequence_str : str
        i : int

        """
        self._sample_aln.replace_sample(sequence_str, i)

    def insert_samples(self, sequence_str, i):
        """Inserts a new sequence in the alignment matrix at the specified
        row position. This increases the total number of rows.

        Parameters
        ----------
        sequence_str : str or list of str
        i : int or list of int

        """
        self._sample_aln.insert_samples(sequence_str, i)

    def append_sample(self, sequence_str):
        """Inserts a new sequence after the last row of the alignment matrix.
        This increases the total number of rows by 1.

        Parameters
        ----------
        sequence_str : str

        """
        self._sample_aln.append_sample(sequence_str)

    def remove_samples(self, i):
        """Removes sample sequences based on the given index.
        If index is a number, only one sequence is removed.
        If the index is a list of numbers, the sequence found at each row
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        self._sample_aln.remove_samples(i)

    def insert_sites(self, sequence_str, i, marker_str=None):
        """Inserts a new sequence in the alignment matrix at the specified
        site position. This increases the total number of columns.

        Parameters
        ----------
        sequence_str : str or list of str
        i : int or list of int

        """
        if marker_str is None and self._marker_aln:
            assert ValueError('marker_str cannot be None if the alignment ' \
                              'has marker sequences')
        if (marker_str is not None) and (not self._marker_aln):
            assert ValueError('The alignment does not use marker sequences')
        self._sample_aln.insert_sites(sequence_str, i)
        if self._marker_aln:
            self._marker_aln.insert_sites(marker_str, i)

    def append_site(self, sequence_str, marker_str=None):
        """Inserts a new sequence at after the last column of the
        alignment matrix. This increases the total number of columns by 1.

        Parameters
        ----------
        sequence_str : str

        """
        if marker_str is None and self._marker_aln:
            assert ValueError('marker_str cannot be None if the alignment ' \
                              'has marker sequences')
        if (marker_str is not None) and (not self._marker_aln):
            assert ValueError('The alignment does not use marker sequences')
        self._sample_aln.append_site(sequence_str)
        if self._marker_aln:
            self._marker_aln.append_site(marker_str)

    def remove_sites(self, i):
        """Removes sites based on the given index.
        If index is a number, only one site is removed.
        If the index is a list of numbers, the sequence found at each column
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        self._sample_aln.remove_sites(i)
        if self._marker_aln:
            self._marker_aln.remove_sites(i)

    def get_samples(self, i):
        """Returns a list of sequence strings containing only the samples
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        return self._sample_aln.get_samples(i)

    def get_markers(self, i):
        """Returns a list of sequence strings containing only the markers
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        return self._marker_aln.get_samples(i)

    def get_sites(self, i):
        """Returns a new alignment containing only the sites specified
        by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        AlignmentMatrix

        """
        return self.__class__.subset(self, sites=i)

    @classmethod
    def from_fasta(cls, path, name=None, marker_kw=None,
                   sample_to_uint_fn=None, uint_to_sample_fn=None,
                   marker_to_uint_fn=None, uint_to_marker_fn=None):
        """Create an Alignment from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file
        marker_kw : str, optional
            A sample is considered a marker if this keyword is present
            within the sequence ID
        sample_to_uint_fn : function, optional
        uint_to_sample_fn : function, optional
        marker_to_uint_fn : function, optional
        uint_to_marker_fn : function, optional

        Raises
        ------
        TypeError

        Returns
        -------
        Alignment

        """
        sequence_list = []
        marker_list = []
        for item in fasta_file_to_list(path, marker_kw=marker_kw):
            if isinstance(item, Sequence):
                sequence_list.append(item)
            elif isinstance(item, Marker):
                marker_list.append(item)
            else:
                raise TypeError('expected Sequence or Marker object')
        return cls(sequence_list, marker_list, name=name,
                   sample_to_uint_fn=sample_to_uint_fn,
                   uint_to_sample_fn=uint_to_sample_fn,
                   marker_to_uint_fn=marker_to_uint_fn,
                   uint_to_marker_fn=uint_to_marker_fn)

    def __repr__(self):
        return '{}(nsamples={}, nsites={}, nmarkers={})'.format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            self.nmarkers
        )

    def __str__(self):
        return '\n'.join([str(self._sample_aln), str(self._marker_aln)])


CatBlock = namedtuple('CatBlock', 'id start stop')


class CatAlignment(Alignment):
    def __init__(self, sequence_list, marker_list, concat_list):
        super().__init__(sequence_list, marker_list)
        self.catblocks = tuple(concat_list)
        self.catblocks_map = {cb.id: self.catblocks[i]
                              for i, cb in enumerate(concat_list)}
        self.block_lists_map = {}

    @classmethod
    def concatenate(cls, aln_list, aln_ids=None, use_aln_names=True):
        start = 0
        def coords(sid, val):
            nonlocal start
            start += val
            return CatBlock(sid, start-val, start)
        # Create a new concat alignment
        new_aln = cls.__new__(cls)
        new_aln.name = 'concat_' + '_'.join([str(aln.name) for aln in aln_list])
        # Create new sample alignment from matrix
        total_sites = sum((aln.nsites for aln in aln_list))
        concat_block_lists = [[Block(0, total_sites)]
                              for i in range(aln_list[0].nsamples)]
        # Put block lists in a mapping
        new_aln.block_lists_map = {
            aln.name: [[Block(b.start, b.stop) for b in blist]
                       for blist in aln.samples.block_lists]
            for aln in aln_list}
        new_aln._sample_aln = SampleAlignment.from_uint_matrix(
            np.concatenate([aln.sample_matrix for aln in aln_list], axis=1),
            aln_list[0].samples.ids,
            [','.join([str(aln.name) for aln in aln_list])] * \
                len(aln_list[0].samples.descriptions),  # replaces desc
            concat_block_lists,  # empties block list
            to_uint_fn=aln_list[0].samples.custom_to_uint_fn,
            from_uint_fn=aln_list[0].samples.custom_from_uint_fn,
            to_block_fn=aln_list[0].samples.custom_to_block_fn,
            from_block_fn=aln_list[0].samples.custom_from_block_fn,
        )
        # Create new marker alignment from matrix
        new_aln._marker_aln = MarkerAlignment.from_uint_matrix(
            np.concatenate([aln.marker_matrix for aln in aln_list], axis=1),
            aln_list[0].markers.ids,
            aln_list[1].markers.descriptions,
            to_uint_fn=aln_list[0].samples.custom_to_uint_fn,
            from_uint_fn=aln_list[0].samples.custom_from_uint_fn,
        )
        if aln_ids is not None:
            new_aln.catblocks = tuple(coords(i, aln.nsites)
                                      for i, aln in zip(aln_ids, aln_list))
        elif use_aln_names:
            new_aln.catblocks = tuple(coords(aln.name, aln.nsites)
                                      for aln in aln_list)
        else:
            new_aln.catblocks = tuple(coords(i, aln.nsites)
                                      for i, aln in enumerate(aln_list))
        new_aln.catblocks_map = {cb.id: new_aln.catblocks[i]
                                 for i, cb in enumerate(new_aln.catblocks)}
        return new_aln

    def get_alignment(self, name):
        if name not in self.catblocks_map:
            raise IndexError("name not found")
        _, start, stop = self.catblocks_map[name]
        aln = Alignment.subset(self, sites=range(start, stop))
        aln.name = 'subaln_{}'.format(name)
        aln.block_lists = [[Block(b.start, b.stop) for b in blist]
                           for blist in self.block_lists_map[name]]
        return aln

    def splitg(self):
        return (Alignment.subset(self, sites=range(cb.start, cb.stop))
                for cb in self.catblocks)

    def split(self):
        return list(self.splitg())

def fasta_file_to_list(path, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str

    Returns
    -------
    list of tuple

    """
    name = ''
    description = ''
    _seq = ''
    seq_list = []
    with open(path, 'r') as f:  # pylint: disable=invalid-name
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                # Store sequence if _seq has contents
                if _seq:
                    if marker_kw:
                        if marker_kw in name:
                            seq = Marker(name, description, _seq)
                        else:
                            seq = Sequence(name, description, _seq)
                    else:
                        seq = Sequence(name, description, _seq)
                    seq_list.append(seq)
                    _seq = ''
                # Split id and description
                try:
                    name, description = line[1:].split(' ', 1)
                except ValueError:
                    name = line[1:]
                    description = ''
            else:
                _seq += line
        if _seq:
            if marker_kw:
                if marker_kw in name:
                    seq = Marker(name, description, _seq)
                else:
                    seq = Sequence(name, description, _seq)
            else:
                seq = Sequence(name, description, _seq)
            seq_list.append(seq)
    return seq_list


def fasta_file_to_alignment(path, marker_kw=None,
                            sample_to_uint_fn=None, uint_to_sample_fn=None,
                            marker_to_uint_fn=None, uint_to_marker_fn=None):
    """Reads a FASTA formatted text file into an Alignment object.

    Parameters
    ----------
    path : str

    Returns
    -------
    list of tuple

    """
    return Alignment.from_fasta(path, marker_kw,
                                sample_to_uint_fn, uint_to_sample_fn,
                                marker_to_uint_fn, uint_to_marker_fn)
