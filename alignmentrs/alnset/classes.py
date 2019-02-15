from collections import OrderedDict
import os
import random
import warnings

from alignmentrs.aln.classes import Alignment, CatAlignment
from libalignmentrs.alignment import concat_basealignments
from libalignmentrs.position import Block, blocks_to_linspace



__all__ = ['AlignmentSet']


class AlignmentMismatchWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass


class DuplicateAlignmentWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass


class AlignmentSet:
    """A container to group Alignments that share the same samples/markers.

    Attributes
    ----------
    name : str
        Name of alignment set.
    metadata : dict
        Metadata about the alignment set.

    """
    def __init__(self, name, aln_list, metadata=None):
        """Creates a new AlignmentSet from a list of Alignment objects.

        Parameters
        ----------
        name : str
            Name of alignment set.
        aln_list : list of Alignment
            List of Alignments to group together.
            Alignments should have thet same number of samples and markers, 
            and share the same sample and marker names.
            A warning is raised when any one of these conditions are
            not satisfied.

        Warnings
        --------
        AlignmentMismatchWarning
            Alignments do not have matching identifiers or have inconsitent markers.

        """
        self.name = name
        self._alignments: OrderedDict = OrderedDict()
        for aln in aln_list:
            if aln.name in self._alignments.keys():
                warnings.warn('Alignment with the same name already exists: {}'
                    .format(aln.name), DuplicateAlignmentWarning)
            self._alignments[aln.name] = aln
        self.metadata: dict = metadata if metadata else dict()
        # Check alignments
        self._consistent: bool = True
        self.drop_empty()


    # Properties
    # ==========================================================================
    @property
    def nalns(self):
        """int: Returns the number of alignments in the set"""
        return len(self)

    @property
    def alignment_names(self):
        """list: Returns list of alignment names"""
        return list(self._alignments.keys())

    @property
    def alignments(self):
        """list of Alignment: Returns list of Alignment objects"""
        return list(self._alignments.values())

    @property
    def consistent(self):
        """bool: Returns whether alignments in the set have consistent
        identifiers and markers"""
        return self._consistent


    # Methods
    # ==========================================================================
    def resample(self, k, with_replacement=False):
        """Takes a sample of alignments from the set and returns a new AlignmentSet. 

        Parameters
        ----------
        k : int
            Number of alignments to sample.
        with_replacement : bool, optional
            Whether to resample with or without replacement.
            If True, resampled alignments may be picked more than once.
            Otherwise, resampled alignments are guaranteed to be unique.
            (default is False)

        Returns
        -------
        list of Alignment
            Note that the alignments are not deep copies of Alignment objects.

        """
        if with_replacement:
            keys = random.choices(self.alignment_names, k=k)
        else:
            keys = random.sample(self.alignment_names, k)
        return [aln for k, aln in self._alignments.items() if k in keys]

    def concatenate(self, name, keys=None):
        """Returns an concatenated alignment from the alignment set.

        Parameters
        ----------
        name : str
            Name of new alignment.
        keys : list or None, optional
            List of alignment names to concatenate. The order of the list
            determines the order of concatenation.
            (default is None, alignments are concatenated arbitrarily)

        Returns
        -------
        Alignment

        """
        if (keys is not None) and not isinstance(keys, list):
            raise ValueError('keys must be None or a list of alignment names')
        # Check if alignments are compatible
        self._check_raise()

        samples = []
        markers = []
        block_list = []
        start = 0
        for k in (name for name in self.alignment_names):
            sample = self._alignments[k].samples
            marker = self._alignments[k].markers
            samples.append(sample)
            markers.append(marker)

            block_list.append(Block(str(k), start, start + sample.nsites))
            start += sample.nsites

        sample_alignment = concat_basealignments(samples)
        if markers[0]:
            marker_alignment = concat_basealignments(markers)
        else:
            marker_alignment = None

        subspaces = OrderedDict({
            str(k): v._linspace for k, v in self._alignments.items()
        })
        return CatAlignment(
            name, sample_alignment, marker_alignment,
            linspace=blocks_to_linspace(block_list),
            subspaces=subspaces,
        )

    def drop_empty(self):
        """Deletes empty alignments and recheck consistency."""
        # Remove empty alignments
        empty_ids = []
        for name, aln in self._alignments.items():
            if not aln:
                empty_ids.append(name)
        for name in empty_ids:
            del self._alignments[name]
        # Check if updated set is consistent
        self._consistent = self._check_warn()


    # Format converters
    # ==========================================================================
    def to_fasta_files(self, path_mapping, include_markers=True,
                       include_headers=True, include_metadata=True):
        """Saves specified alignments as FASTA-formatted text files
        based on a dictionary of alignment names and their corresponding
        save paths.

        Parameters
        ----------
        path_mapping : dict
            Keys are alignment names, values are paths where to
            write alignments.
        include_markers : bool, optional
            If True, markers are also written to the FASTA file.
            (default is True)
        include_headers : bool, optional
            Whether or not to output header infomation,
            ie. alignment name and coordinates.
            (default is True, include headers in the output as
            comments)
        include_metadata : bool, optional
            Whether or not to output metadata as comments.
            (default is True, include metadata in the output as
            comments)

        Returns
        -------
        int
            Number of FASTA files written.

        """
        c = 0
        for key, path in path_mapping.items():
            if key not in self._alignments.keys():
                raise KeyError('path_mapping key "{}" does not match '
                               'any alignment name.'.format(key))
            self._alignments[key].to_fasta(
                path,
                include_markers=include_markers,
                include_headers=include_headers,
                include_metadata=include_metadata)
            c += 1
        return c

    def to_fasta_dir(self, dirpath, include_markers=True,
                     include_headers=True, include_metadata=True,
                     name_to_filename_encoder=None):
        """Saves all alignments as FASTA-formatted text files
        in a single directory.

        Parameters
        ----------
        dirpath : str
            Folder/directory to save the FASTA files in.
        include_markers : bool, optional
            If True, markers are also written to the FASTA file.
            (default is True)
        include_headers : bool, optional
            Whether or not to output header infomation,
            ie. alignment name and coordinates.
            (default is True, include headers in the output as
            comments)
        include_metadata : bool, optional
            Whether or not to output metadata as comments.
            (default is True, include metadata in the output as
            comments)
        name_to_filename_encoder : function or None, optional
            This function returns the alignment's filename
            given its name. The function receives one parater,
            the alignment name.
            If not specified, the output filename will be in the
            format `{name}.aln` where name is the alignment's name.
            (default is None)

        Returns
        -------
        int
            Number of FASTA files written.

        """
        if not os.path.exists(dirpath):
            os.makedirs(os.path.abspath(dirpath))
        if name_to_filename_encoder is None:
            name_to_filename_encoder = lambda x: '{}.aln'.format(x)
        path_mapping = {
            name: os.path.join(dirpath, name_to_filename_encoder(name))
            for name in self._alignments.keys()}
        return self.to_fasta_files(
            path_mapping,
            include_markers=include_markers,
            include_headers=include_headers,
            include_metadata=include_metadata)

    @classmethod
    def from_fasta_files(cls, paths, name, marker_kw=None,
                         filename_to_key_encoder=None,
                         comment_parser=None):
        """Reads FASTA files and stores data as a set of Alignment objects
        inside AlignmentSet.

        Parameters
        ----------
        paths : list
            List of FASTA file paths.
        name : str
            Name of the alignment set.
        marker_kw : str or None, optional
            Classifies the sample as a marker if the string is
            found in the sample's ID. (default is None)
        filename_to_key_encoder : function or None, optional
            If specified, the function receives the filename as input
            and outputs a key to identify a unique alignment.
            This can be used to make sure that the same alignment
            stored as files with different filenames are not
            included multiple times.

        Returns
        -------
        AlignmentSet
            New AlignmentSet object with each FASTA file as a member Alignment
            object.

        """
        sequence_d = {}
        for fname in paths:
            key = filename_to_key_encoder(fname) \
                  if filename_to_key_encoder else fname
            if key in sequence_d.keys():
                raise KeyError('alignment "{}" already exists'.format(key))
            sequence_d[key] = Alignment.from_fasta(
                fname, name=key, marker_kw=marker_kw,
                comment_parser=comment_parser)
        return cls(name, list(sequence_d.values()))

    @classmethod
    def from_fasta_dir(cls, dirpath, name, marker_kw=None, suffix='.aln',
                       filename_to_key_encoder=None,
                       comment_parser=None):
        """Reads a directory containing FASTA files and stores data as a
        set of alignment objects inside an AlignmentSet.

        Parameters
        ----------
        dirpath : str
            Path containing FASTA files to be read.
        name : str
            Name of alignment set.
        marker_kw : str or None, optional
            Classifies the sample as a marker if the string is
            found in the sample's ID. (default is None)
        suffix : str, optional
            Used to determine whether a file is a FASTA file (default is '.aln')
        marker_kw : str or None, optional
            Classifies the sample as a marker if the string is
            found in the sample's ID. (default is None)
        filename_to_key_encoder : function or None, optional
            If specified, the function receives the filename as input
            and outputs a key to identify a unique alignment.
            This can be used to make sure that the same alignment
            stored as files with different filenames are not
            included multiple times.

        Returns
        -------
        AlignmentSet
            New AlignmentSet object with each FASTA file as a member Alignment
            object.

        """
        # Check if dirpath exists
        if not os.path.exists(dirpath):
            raise Exception('{} does not exist'.format(dirpath))
        else:
            if not os.path.isdir(dirpath):
                raise Exception('{} is not a directory'.format(dirpath))
        paths = (os.path.join(dirpath, fname) for fname in os.listdir(dirpath)
                 if fname.endswith(suffix))
        return cls.from_fasta_files(
            paths, name, marker_kw,
            filename_to_key_encoder,
            comment_parser)


    # Check alignments
    # ==========================================================================
    def _check_warn(self):
        # Check if alignments are compatible
        sample_missing = 0
        sample_mismatch = 0

        marker_present = 0
        marker_missing = 0
        marker_mismatch = 0

        test_aln = None        
        for _, aln in self._alignments.items():
            if not aln.samples:
                sample_missing += 1
            if aln.markers:
                marker_present += 1
            else:
                marker_missing += 1

            if test_aln is None:
                test_aln = aln
                continue

            if test_aln.samples.nrows != aln.samples.nrows:
                sample_mismatch += 1
            if test_aln.markers.nrows != aln.markers.nrows:
                marker_mismatch += 1

            test_aln = aln

        passed = True
        if sample_missing > 0:
            msg = '{}/{} alignments do not have sample entries.'.format(
                sample_missing, len(self._alignments)
            )
            warnings.warn(msg, AlignmentMismatchWarning)
            passed = False

        if sample_mismatch > 0:
            msg = '{}/{} alignments do not have the same number ' \
                  'of samples.'.format(sample_mismatch, len(self._alignments)-1)
            warnings.warn(msg, AlignmentMismatchWarning)
            passed = False

        if marker_present > 0 and marker_missing > 0:
            msg = '{}/{} alignments have markers, {}/{} alignments do not have markers.'.format(
                marker_present, len(self._alignments), 
                marker_mismatch, len(self._alignments), 
            )
            warnings.warn(msg, AlignmentMismatchWarning)
            passed = False

        if marker_mismatch > 0:
            msg = '{}/{} alignments do not have the same number ' \
                  'of markers.'.format(marker_mismatch, len(self._alignments)-1)
            warnings.warn(msg, AlignmentMismatchWarning)
            passed = False

        return passed

    def _check_raise(self):
        test_aln = None
        for _, aln in self._alignments.items():
            if test_aln is None:
                test_aln = aln
                continue

            if not(test_aln.samples or aln.samples):
                raise ValueError('Alignment is missing sample sequences.')
            else:
                if test_aln.samples.nrows != aln.samples.nrows:
                    raise ValueError('Number of samples do not match.')

            if not(test_aln.markers != aln.markers):
                raise ValueError('Presence/absence of markers is not consistent.')
            else:
                if test_aln.markers.nrows != aln.markers.nrows:
                    raise ValueError('Number of markers do not match')


    # Special methods
    # ==========================================================================
    def __getitem__(self, key):
        return self._alignments[key]

    def __delitem__(self, key):
        self._alignments.__delitem__(key)

    def __iter__(self):
        for name, aln in self._alignments.items():
            yield name, aln

    def __repr__(self):
        nsamples = 'Inconsistent'
        nmarkers = 'Inconsistent'
        if self.consistent:
            nsamples = self.alignments[0].nsamples if self.nalns > 0 else 'None'
            nmarkers = self.alignments[0].nmarkers if self.nalns > 0 else 'None'
        return '{}(nalns={}, nsamples={}, nmarkers={})'.format(
            self.__class__.__name__, self.nalns, nsamples, nmarkers
        )

    def __len__(self):
        return len(self._alignments)
