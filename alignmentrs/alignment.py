import itertools
import numpy as np
from alignmentrs import Sequence, Marker


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
        self.to_uint_fn = np.vectorize(ord)
        self.from_uint_fn = np.vectorize(chr)
        if to_uint_fn is not None:
            self.to_uint_fn = np.vectorize(to_uint_fn)
        if from_uint_fn is not None:
            self.from_uint_fn = np.vectorize(from_uint_fn)
        self.matrix = np.empty((nsamples, nsites), dtype=np.uint32)

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
        mat = cls.__new__(cls)
        mat.matrix = np.empty((nsamples, nsites), dtype=np.uint32)
        return mat

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
        new_aln_matrix = cls(*new_matrix.shape)
        new_aln_matrix.matrix = new_matrix
        return new_aln_matrix

    @classmethod
    def from_matrix(cls, matrix, to_uint_fn=None, from_uint_fn=None):
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

        new_matrix = cls(*matrix.shape,
                         to_uint_fn=to_uint_fn,
                         from_uint_fn=from_uint_fn)
        for i, row in enumerate(matrix):
            new_matrix.replace_sample(''.join(row), i)
        
        return new_matrix

    def replace_sample(self, i, sequence_str):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        sequence_str : str
        i : int

        """
        self.matrix[i, :] = [self.to_uint_fn(s) for s in sequence_str]

    def insert_sample(self, i, sequence_str):
        """Inserts a new sequence in the alignment matrix at the specified
        row position. This increases the total number of rows.

        Parameters
        ----------
        sequence_str : str
        i : int

        """
        new_array = np.array([[self.to_uint_fn(s) for s in sequence_str]],
                             dtype=np.uint32)
        self.matrix = np.insert(self.matrix, i, new_array, axis=0)

    def append_sample(self, sequence_str):
        """Inserts a new sequence after the last row of the alignment matrix.
        This increases the total number of rows by 1.

        Parameters
        ----------
        sequence_str : str

        """
        new_array = np.array([[self.to_uint_fn(s) for s in sequence_str]],
                             dtype=np.uint32)
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
        translated_g = (self.from_uint_fn(row) for row in self.matrix[i, :])
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
        translated_g = (self.from_uint_fn(row) for row in self.matrix[:, i])
        return [''.join(trans) for trans in translated_g]

    def replace_site(self, i, sequence_str):
        """Replaces the sequence for a given column in the alignment matrix.

        Parameters
        ----------
        sequence_str : str
        i : int

        """
        self.matrix[:, i] = [self.to_uint_fn(s) for s in sequence_str]

    def insert_site(self, i, sequence_str):
        """Inserts a new sequence in the alignment matrix at the specified
        site position. This increases the total number of columns.

        Parameters
        ----------
        sequence_str : str
        i : int

        """
        new_array = np.array([self.to_uint_fn(s) for s in sequence_str],
                             dtype=np.uint32)
        self.matrix = np.insert(self.matrix, i, new_array, axis=1)

    def append_site(self, sequence_str):
        """Inserts a new sequence at after the last column of the
        alignment matrix. This increases the total number of columns by 1.

        Parameters
        ----------
        sequence_str : str

        """
        new_array = np.array([[self.to_uint_fn(s) for s in sequence_str]],
                             dtype=np.uint32).T
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

    def __len__(self):
        return len(self.matrix)

    def __repr__(self):
        return "{}(nsamples={}, nsites={})\n{}".format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            repr(self.matrix)
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
        super().__init__(len(sequence_list), len(sequence_list[0].sequence),
                         to_uint_fn=to_uint_fn, from_uint_fn=from_uint_fn)
        self.ids = [s.id for s in sequence_list]
        self.descriptions = [s.description for s in sequence_list]
        # TODO: Adds option to generate different matrix if
        # to_uint_fn, from_uint_fn is not None
        self.matrix = np.array([s.sequence_to_uint32() for s in sequence_list])

    @classmethod
    def subset(cls, aln, rows=None, cols=None,
               row_step=1, col_step=1):
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
        new_aln = cls.__new__(cls)
        new_aln.ids = [sid for i, sid in enumerate(aln.ids) if i in rows]
        new_aln.descriptions = [desc for i, desc in enumerate(aln.descriptions)
                                if i in rows]
        new_aln.matrix = np.copy(aln.matrix[rows][:, cols])

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
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            elif len(i) == sum((isinstance(j, str) for j in i)):
                # get int position in self.ids and then use subset
                i = [self.ids.index(idx) for idx in i]
            else:
                raise ValueError('i must be a list of int or str')
        else:
            raise ValueError('i must be an int or a list of int or str')
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
        translated_g = (self.from_uint_fn(row) for row in sample_array.matrix)
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
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            elif len(i) == sum((isinstance(j, str) for j in i)):
                # get int position in self.ids and then use subset
                i = [self.ids.index(idx) for idx in i]
            else:
                raise ValueError('i must be a list of int or str')
        else:
            raise ValueError('i must be an int or a list of int or str')
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
        elif isinstance(i, list):
            if len(i) == sum((isinstance(j, int) for j in i)):
                return self.__class__.subset(self, rows=i)
            else:
                raise ValueError('i must be a list of int')
        else:
            raise ValueError('i must be an int or a list of int')
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
        translated_g = (self.from_uint_fn(row) for row in sample_array.matrix)
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

    def __str__(self):
        id_lines = ['>{} {}'.format(sid, desc) if desc else '>{}'.format(sid)
                    for sid, desc in zip(self.ids, self.descriptions)]
        seq_lines = self.get_samples_as_str(list(range(len(self))))
        return '\n'.join(itertools.chain(*zip(id_lines, seq_lines)))


class SampleAlignment(BaseAlignment):
    """SampleAlignment represents a mulitple sequence alignment of
    biological sequences.
    """
    pass


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
        super().__init__(marker_list, to_uint_fn=to_uint_fn, from_uint_fn=from_uint_fn)

class Alignment(object):
    """Alignment is a complete representation of a multiple sequence alignment
    of biological sequences an their annotations such as alignment markers and
    alignment block data.
    """
    def __init__(self, sequence_list, marker_list,
                 sequence_to_uint_fn=None, uint_to_sequence_fn=None,
                 marker_to_uint_fn=None, uint_to_marker_fn=None):
        """Creates a new Alignment object from a list of Sequence and Marker objects.

        Parameters
        ----------
        sequence_list : list of Sequence
        marker_list : list of Marker

        """
        self._sample_aln = SampleAlignment(sequence_list,
                                           sequence_to_uint_fn,
                                           uint_to_sequence_fn)
        self._marker_aln = MarkerAlignment(marker_list,
                                           marker_to_uint_fn,
                                           uint_to_marker_fn)
        assert self._sample_aln.nsites == self._marker_aln.nsites

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
    def from_fasta(cls, path, marker_kw='',
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

        sequences = [s for s in fasta_file_to_list(path)]
        item_list = fasta_file_to_list(path, marker_kw=marker_kw)
        sequence_list = []
        marker_list = []
        for item in item_list:
            if isinstance(item, Sequence):
                sequence_list.append(item)
            elif isinstance(item, Marker):
                marker_list.append(item)
            else:
                raise TypeError('expected Sequence or Marker object')
        return cls(sequence_list, marker_list,
                   sample_to_uint_fn, uint_to_sample_fn,
                   marker_to_uint_fn, uint_to_marker_fn)

    def __repr__(self):
        return '{}(nsamples={}, nsites={}, nmarkers={})'.format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            self.nmarkers
        )

    def __str__(self):
        return '\n'.join([str(self._sample_aln), str(self._marker_aln)])

class BlockAlignment(Alignment):
    def __init__(self, sequence_list, marker_list, block_data_parser=None):
        super().__init__(sequence_list, marker_list)
        self.custom_block_data_parser = None

    @staticmethod
    def parse_block_data(string):
        pass

class ConcatAlignment(Alignment):
    def __init__(self, sequence_list, marker_list, concat_data_parser=None):
        super().__init__(sequence_list, marker_list)
        self.concat_data_parser = self.parse_concat_data \
            if concat_data_parser is None else concat_data_parser

    @staticmethod
    def parse_concat_data(string):
        pass


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
            if marker_kw is not None:
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
