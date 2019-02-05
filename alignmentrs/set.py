from alignmentrs.alignment import Alignment
from libalignmentrs.alignment import concat_basealignments
import random
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
