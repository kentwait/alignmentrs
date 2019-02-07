import os
import random

from blockrs import CatBlock
from libalignmentrs.alignment import concat_basealignments
from alignmentrs.aln.classes import Alignment


__all__ = ['AlignmentSet']


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
            List of Alignments to group together. Alignments must have thet
            same number of samples and markers, and share the same sample
            and marker names.

        Raises
        ------
        ValueError
            Alignments in `aln_list` do not have the same number of
            samples/markers, or so not share the same sample/marker names.

        """
        # Check
        for i in range(1, len(aln_list)):
            if aln_list[0].samples and aln_list[i].samples:
                pass
            elif aln_list[0].markers and aln_list[i].markers:
                pass
            elif aln_list[0].samples.is_row_similar(aln_list[i].samples):
                pass
            elif aln_list[0].markers.is_row_similar(aln_list[i].markers):
                pass
            else:
                raise ValueError(
                    'Cannot create an AlignmentSet from a alignment with '
                    'different number of samples/markers, or different '
                    'sample/marker names.')
        self.name = name
        self._alignments = {aln.name: aln for aln in aln_list}
        self.metadata = metadata if metadata else dict()

    @property
    def nalns(self):
        """int: Returns the number of alignments in the set"""
        return len(self)

    @property
    def alignment_names(self):
        """list of str: Returns list of alignment names"""
        return list(self._alignments.keys())

    @property
    def alignments(self):
        """list of Alignment: Returns list of alignment objects"""
        return list(self._alignments.values())

    def resample(self, name, k, copy_metadata=True,
                 with_replacement=False):
        """Returns a new AlignmentSet resampling from the current
        set of alignments.

        Parameters
        ----------
        name : str
            Name of resampled alignment set.
        k : int
            Number of alignments to sample.
        copy_metadata : bool, optional
            Copies existing metadata to the resampled set if True.
        with_replacement : bool, optional
            Whether to resample with or without replacement.
            If True, resampled alignments may be picked more than once.
            Otherwise, resampled alignments are guaranteed to be unique.
            (default is False)

        Returns
        -------
        AlignmentSet
            Creates a new AlignmentSet grouping. Note that the alignments
            are not deep copies as AlignmentSet a grouping construct only.

        Notes
        -----
        The resampled alignment set will have an additional metadata entry
        with the key "resampling", if it is not present already.
        This entry will have a dictionary containing the following
        resampling information:
        - number of original items as `nalns`
        - number of sampled alignments as `k`
        - with/without replacement as `with_replacement`

        If resampling already exists, it will be overwritten by the new
        values.

        """
        if with_replacement:
            keys = random.choices(self.alignment_names, k=k)
        else:
            keys = random.sample(self.alignment_names, k)
        aln_list = [aln for k, aln in self._alignments.items() if k in keys]
        if copy_metadata:
            new_metadata = {k: v for k, v in self.metadata.items()}
        else:
            new_metadata = dict()
        new_metadata['resampling'] = {
            'nalns': len(self._alignments),
            'k': k,
            'with_replacement': with_replacement,
        }
        return self.__class__(name, aln_list, metadata=new_metadata)

    def concatenate(self, name, keys=None, description_encoder=None):
        """Returns an concatenated alignment from the alignment set.

        Parameters
        ----------
        name : str
            Name of new alignment.
        keys : list or None, optional
            List of alignment names to concatenate. The order of the list
            determines the order of concatenation.
            If None, alignments are concatenated based on original order.
            (default is None)
        description_encoder : function or None, optional
            This function returns a formatted string encoding CatBlock data.
            The CatBlock string will replace the sample's description.
            This function receives two parameters, the sample ID and the
            sample's concatenation order as a list of CatBlocks.
            If None, the description will follow the description of the
            first alignment in during concatenation.

        Returns
        -------
        Alignment

        """
        if keys is not None:
            if not isinstance(keys, list):
                raise ValueError('keys must be None or a list of alignment names')
        else:
            keys = (name for name in self.alignment_names)
        sample_alignment = concat_basealignments(
            [self._alignments[k].samples for k in keys])
        try:
            marker_alignment = concat_basealignments(
                [self._alignments[k].markers for k in keys])
        except Exception:
            marker_alignment = None

        if description_encoder:
            start = 0
            def cb(sid, val):
                nonlocal start
                start += val
                return CatBlock(str(sid), start-val, start)
            cblist = [cb(name, aln.nsites)
                      for name, aln in self._alignments.items()]
            descriptions = [description_encoder(sid, cblist)
                            for sid in sample_alignment.ids]
            sample_alignment.set_descriptions(
                list(range(sample_alignment.nsamples)),
                descriptions
            )
        return Alignment(name, sample_alignment, marker_alignment)

    def to_fasta_files(self, path_mapping, include_markers=True):
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
            self._alignments[key].to_fasta(path, include_markers)
            c += 1
        return c

    def to_fasta_dir(self, dirpath, include_markers=True,
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
        return self.to_fasta_files(path_mapping, include_markers)

    @classmethod
    def from_fasta_files(cls, paths, name, marker_kw=None,
                         filename_to_key_encoder=None):
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
            sequence_d[key] = Alignment.from_fasta(fname, name=key,
                                                   marker_kw=marker_kw)
        return cls(name, list(sequence_d.values()))

    @classmethod
    def from_fasta_dir(cls, dirpath, name, marker_kw=None,
                       suffix='.aln', filename_to_key_encoder=None):
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
        return cls.from_fasta_files(paths, name, marker_kw,
                                    filename_to_key_encoder)

    def __getitem__(self, key):
        return self._alignments[key]

    def __delitem__(self, key):
        self._alignments.__delitem__(key)

    def __iter__(self):
        for name, aln in self._alignments.items():
            yield name, aln

    def __repr__(self):
        return '{}(nalns={}, nsamples={}, nmarkers={})'.format(
            self.__class__.__name__,
            self.nalns,
            self.alignments[0].nsamples if self.nalns > 0 else None,
            self.alignments[0].nmarkers if self.nalns > 0 else None
        )

    def __len__(self):
        return len(self._alignments)
