from alignmentrs.alignment import Alignment
from libalignmentrs.alignment import concat_basealignments
import random
import os
import blockrs


class AlignmentSet:
    def __init__(self, name, aln_list):
        self.name = name
        self._alignments = {aln.name: aln for aln in aln_list}

    @property
    def nalns(self):
        """Returns the number of alignments in the set"""
        return len(self)

    @property
    def alignment_names(self):
        """Returns list of alignment names"""
        return list(self._alignments.keys())

    @property
    def alignments(self):
        """Returns list of alignment objects"""
        return list(self._alignments.values())

    def resample(self, k, name='resampled_alignment_set',
                 with_replacement=False):
        """Returns a new AlignmentSet resampling from the current
        set of alignments.

        Parameters
        ----------
        k : int
            Number of alignments
        name : str, optional
            Name of resampled superalignment. (default is 'resampled_alignment_set')
        with_replacement : bool, optional
            If True, resampled alignments may not be unique.
            If False, resampled alignments are guaranteed to be unique.

        """
        if with_replacement:
            keys = random.choices(self.alignment_names, k=k)
        else:
            keys = random.sample(self.alignment_names, k)
        aln_list = [aln for k, aln in self._alignments.items() if k in keys]
        return self.__class__(name, aln_list)

    def concatenate(self, name='concatenated_alignment', keys=None):
        """Returns an concatenated alignment from the alignment set.

        Parameters
        ----------
        name : str, optional
            Name of new alignment. (default is 'concatenated_alignment')
        keys : list or None, optional
            List of alignment names to concatenate. The order of the list
            determines the order of concatenation.
            If None, alignments are concatenated based on original order.
            (default is None)

        Returns
        -------
        Alignment

        """
        if keys is not None:
            if isinstance(keys, list) and \
               sum((isinstance(k, str) for k in keys)):
                pass
            else:
                raise ValueError('keys must be None or a list of alignment names')
        sample_alignment = concat_basealignments(
            [aln.samples for aln in self._alignments.values()])
        try:
            marker_alignment = concat_basealignments(
                [aln.markers for aln in self._alignments.values()])
        except Exception:
            marker_alignment = None
        return Alignment(name, sample_alignment, marker_alignment)

    @classmethod
    def from_fasta_files(cls, paths, name, marker_kw=None,
                         filename_to_key_encoder=None):
        """Reads FASTA files and stores data as a set of Alignment objects
        inside AlignmentSet.

        Parameters
        ----------
        paths : list
            List of FASTA file paths
        name : str
            Name of the alignment set
        marker_kw : str or None, optional
            Classifies the sample as a marker if the string is
            found in the sample's ID. (default is None)
        filename_to_key_encoder : function or None, optional
            If specified, the function receives the filename as input
            and outputs a key to identify a unique alignment.
            THis can be used to make sure that the same alignment
            stored as files with different filenames are not
            included multiple times.

        Returns
        -------
        AlignmentSet

        """
        sequence_d = {}
        for fname in paths:
            key = filename_to_key_encoder(fname) \
                  if filename_to_key_encoder else fname
            if key in sequence_d.keys():
                raise KeyError('alignment "{}" already exists'.format(key))
            sequence_d[key] = Alignment.from_fasta(fname, name=key,
                                                   marker_kw=marker_kw)
        return cls(name, sequence_d.values())

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
            Name of alignment set
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
            THis can be used to make sure that the same alignment
            stored as files with different filenames are not
            included multiple times.

        Returns
        -------
        AlignmentSet

        """
        # Check if dirpath exists
        if not os.path.exists(dirpath):
            raise Exception('{} does not exist'.format(dirpath))
        else:
            if not os.path.isdir(dirpath):
                raise Exception('{} is not a directory'.format(dirpath))
        paths = (fname for fname in os.listdir(dirpath)
                 if fname.endswith(suffix))
        return cls.from_fasta_files(paths, name, marker_kw,
                                    filename_to_key_encoder)

    def __getitem__(self, key):
        return self._alignments[key]

    def __delitem__(self, key):
        del self._alignments[key]

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


def fasta_directory_to_alignmentset(dirpath, name, marker_kw=None,
                                    suffix='.aln',
                                    filename_to_key_encoder=None):
    """Reads a directory containing FASTA files and stores data as a
    set of alignment objects inside an AlignmentSet.

    Parameters
    ----------
    dirpath : str
        Path containing FASTA files to be read.
    name : str
        Name of alignment set
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
        THis can be used to make sure that the same alignment
        stored as files with different filenames are not
        included multiple times.

    Returns
    -------
    AlignmentSet

    """
    return AlignmentSet.from_fasta_dir(dirpath, name, marker_kw, suffix,
                                       filename_to_key_encoder)
