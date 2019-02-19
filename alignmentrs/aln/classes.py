from collections import OrderedDict
import itertools
import os
from copy import deepcopy

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.readers import (fasta_file_to_basealignment, 
                                    fasta_file_to_basealignments)
from libalignmentrs.position import BlockSpace
from libalignmentrs.record import Record
from alignmentrs.util import parse_comment_list, parse_cat_comment_list
from .mixins import PropsMixin, AlnMixin
from .mixins import SamplePropsMixin, SampleAlnMixin
from .mixins import MarkerPropsMixin, MarkerAlnMixin
from .mixins import FastaSerde, JsonSerde
from .mixins import CoordsMixin
from .mixins.functions import subset


__all__ = ['Alignment', 'CatAlignment']


class Alignment(CoordsMixin, AlnMixin, PropsMixin, object):
    """Reperesents a multiple sequence alignment of samples.

    The Alignment object encapsulates information generally
    included in the FASTA format:
    - sequence names/ids
    - descriptions
    - sequences

    Additionally, the Alignment object also stores comments
    (lines prefixed by a semicolon ";") as metadata.

    Attributes
    ----------
    name : str
        Name of the alignment.
    samples : BaseAlignment
        Alignment of sample sequences.
    metadata : dict
        Other information related to the alignment.
    """

    members = ['samples']

    def __init__(self, name, sample_alignment: BaseAlignment, 
                 linspace: BlockSpace=None, metadata: dict=None,
                 **kwargs):
        """Creates a new Alignment object from a sample BaseAlignment.

        Parameters
        ----------
        name : str
            Name of the alignment.
        sample_alignment : BaseAlignment
            Alignment of sample sequences.
        linspace : BlockSpace, optional
            Linear space that assigns coordinate values to alignment
            columns. (default is None, the constructor will create a new
            linear space starting from 0).
        metadata : dict, optional
            Other information related to the alignment. (default is None,
            which creates a blank dictionary)
        **kwargs
            Other keyword arguments used to initialize states in the
            Alignment object.

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
        self.name = name
        self.samples: BaseAlignment = sample_alignment

        # Use given metadata
        if metadata is not None:
            if isinstance(metadata, dict):
                self.metadata: OrderedDict = OrderedDict(metadata)
            else:
                raise TypeError('metadata is not a dictionary.')
        else:
            self.metadata: OrderedDict = OrderedDict()

        # Use given linspace
        if linspace is not None:
            if isinstance(linspace, BlockSpace):
                self._linspace: BlockSpace = linspace
            else:
                raise TypeError('linspace is not a BlockSpace object.')
        else:
            start = kwargs['linspace_default_start'] \
                    if 'linspace_default_start' in kwargs.keys() else 0
            stop = kwargs['linspace_default_stop'] \
                    if 'linspace_default_stop' in kwargs.keys() else \
                    start + self.samples.ncols
            state = kwargs['linspace_default_state'] \
                    if 'linspace_default_state' in kwargs.keys() else \
                    "1"
            self._linspace: BlockSpace = BlockSpace(start, stop, state)


    @classmethod
    def _from_basealignment(cls, name: str, sample_alignment: BaseAlignment,
                            linspace: BlockSpace=None, metadata: dict=None,
                            **kwargs):
        obj = super(Alignment, cls).__new__(cls)
        if not sample_alignment:
            raise ValueError(
                'Cannot create an Alignment using an empty '
                'sample_alignment BaseAlignment.')
        obj.name = name
        obj.samples: BaseAlignment = sample_alignment

        # Use given metadata
        if metadata is not None:
            if isinstance(metadata, dict):
                obj.metadata: OrderedDict = OrderedDict(metadata)
            else:
                raise TypeError('metadata is not a dictionary.')
        else:
            obj.metadata: OrderedDict = OrderedDict()

        # Use given linspace
        if linspace is not None:
            if isinstance(linspace, BlockSpace):
                obj._linspace: BlockSpace = linspace
            else:
                raise TypeError('linspace is not a BlockSpace object.')
        else:
            start = kwargs['linspace_default_start'] \
                    if 'linspace_default_start' in kwargs.keys() else 0
            stop = kwargs['linspace_default_stop'] \
                    if 'linspace_default_stop' in kwargs.keys() else \
                    start + obj.samples.ncols
            state = kwargs['linspace_default_state'] \
                    if 'linspace_default_state' in kwargs.keys() else \
                    "1"
            obj._linspace: BlockSpace = BlockSpace(start, stop, state)
        return obj


    # Special methods
    # ==========================================================================
    def copy(self):
        """Creates a new deep copy of the alignment.

        Returns
        -------
        Alignment
            The new alignment object is a deep copy of the original alignment.

        """
        new_alns = (self.__getattribute__(member).copy()
                    for member in self.__class__.members)
        return self.__class__._from_basealignment(
            self.name, *new_alns,
            metadata=deepcopy(self.metadata),
            linspace=self._linspace.copy())

    def __getitem__(self, key):
        if isinstance(key, str):
            for member in self.__class__.members:
                if key in self.__getattribute__(member).ids:
                    i = self.__getattribute__(member).row_names_to_ids([key])
                    assert len(i) == 1, \
                        '{} matched multiple rows in the alignment.'.format(key)
                    return self.samples.get_row(i[0])
            raise KeyError('Key did not match any sample name/identifier')
        elif isinstance(key, int):  # TODO: Fix bug
            return self.samples.get_col(key)
        raise TypeError('Key must be str or int.')

    def __delitem__(self, key):
        if isinstance(key, str):
            for member in self.__class__.members:
                if key in self.__getattribute__(member).ids:
                    i = self.samples.row_names_to_ids([key])
                    assert len(i) == 1, \
                        '{} matched multiple rows in the alignment.'.format(key)
                    return self.__getattribute__(member).remove_rows(i)
            raise KeyError('Key did not match any sample name/identifier')
        elif isinstance(key, int):
            return self.remove_cols(key)
        raise TypeError('Key must be str or int.')

    def __iter__(self):
        raise NotImplementedError(
            'Iteration over alignment is ambiguous.\n'
            'Use .iter_samples, .iter_markers, or .iter_rows '
            'to iterate across samples, markers or all entries, respectively.\n'
            'Use .iter_sample_sites, .iter_marker_sites, .iter_sites '
            'to iterate across sites in samples, markers or '
            'sites across both samples and markers, respectively.')

    def __repr__(self):
        if len(self.__class__.members) > 1:
            member_reprs = ', '.join([
                'n{}={}'.format(member, self.__getattribute__(member).nrows) 
                for member in self.__class__.members
            ])
        else:
            member_reprs = 'nrows={}'.format(self.nrows)
        return '{}({}, ncols={})'.format(
            self.__class__.__name__,
            member_reprs,
            self.ncols,
        )

    def __str__(self):
        return self.to_fasta_str()

    def __len__(self):
        raise NotImplementedError(
            'len() is not implemented for Alignment.\n'
            'Use .ncols to get the number of columns, '
            '.nsamples to get the number of samples, '
            '.nmarkers to get the number of markers, or '
            '.nrows to get all the number of alignment rows.')

    def __bool__(self):
        if self.ncols == 0 or self.nrows == 0:
            return False
        return True

    def __hash__(self):
        return hash((
            str(self),
            self._linspace.to_block_str(),
            tuple(self.metadata.items())
        ))

    def __eq__(self, other):
        return hash(self) == hash(other)


class SampleAlignment(JsonSerde, FastaSerde, SampleAlnMixin, SamplePropsMixin, Alignment):
    members = ['samples']


class MarkerAlignment(JsonSerde, FastaSerde, MarkerAlnMixin, MarkerPropsMixin, Alignment):
    members = ['markers']

    @classmethod
    def _from_basealignment(cls, name, marker_alignment: BaseAlignment, 
                            linspace: BlockSpace=None, metadata: dict=None,
                            **kwargs):
        obj = super(MarkerAlignment, cls).__new__(cls)
        if not marker_alignment:
            raise ValueError(
                'Cannot create an Alignment using an empty '
                'marker_alignment BaseAlignment.')
        obj.name = name
        obj.markers: BaseAlignment = marker_alignment

        # Use given metadata
        if metadata is not None:
            if isinstance(metadata, dict):
                obj.metadata: OrderedDict = OrderedDict(metadata)
            else:
                raise TypeError('metadata is not a dictionary.')
        else:
            obj.metadata: OrderedDict = OrderedDict()

        # Use given linspace
        if linspace is not None:
            if isinstance(linspace, BlockSpace):
                obj._linspace: BlockSpace = linspace
            else:
                raise TypeError('linspace is not a BlockSpace object.')
        else:
            start = kwargs['linspace_default_start'] \
                    if 'linspace_default_start' in kwargs.keys() else 0
            stop = kwargs['linspace_default_stop'] \
                    if 'linspace_default_stop' in kwargs.keys() else \
                    start + obj.samples.ncols
            state = kwargs['linspace_default_state'] \
                    if 'linspace_default_state' in kwargs.keys() else \
                    "1"
            obj._linspace: BlockSpace = BlockSpace(start, stop, state)
        return obj


class FullAlignment(JsonSerde, FastaSerde, MarkerAlnMixin, MarkerPropsMixin,
                    SampleAlnMixin, SamplePropsMixin, Alignment):
    """Represents a multiple sequence alignment.

    The Alignment object encapsulates information generally
    included in the FASTA format:
    - sequence names/ids
    - descriptions
    - sequences

    Additionally, the Alignment object also stores comments
    (lines prefixed by a semicolon ";") as metadata.

    Attributes
    ----------
    name : str
        Name of the alignment.
    samples : BaseAlignment
        Alignment of sample sequences.
    markers : BaseAlignment
        Alignment of non-sample sequences. This is metadata
        stored as a row in the alignment that describes some
        kind of site-specific information.
    metadata : dict
        Other information related to the alignment.

    """
    members = ['samples', 'markers']

    @classmethod
    def _from_basealignment(cls, name, 
                            sample_alignment: BaseAlignment,
                            marker_alignment: BaseAlignment, 
                            linspace: BlockSpace=None, metadata: dict=None,
                            **kwargs):
        obj = super(MarkerAlignment, cls).__new__(cls)        
        if not sample_alignment:
            raise ValueError(
                'Cannot create a FullAlignment object using an empty '
                'sample_alignment BaseAlignment.')
        obj.name = name
        obj.samples: BaseAlignment = sample_alignment
        obj.markers: BaseAlignment = marker_alignment

        # Use given metadata
        if metadata is not None:
            if isinstance(metadata, dict):
                obj.metadata: OrderedDict = OrderedDict(metadata)
            else:
                raise TypeError('metadata is not a dictionary.')
        else:
            obj.metadata: OrderedDict = OrderedDict()

        # Use given linspace
        if linspace is not None:
            if isinstance(linspace, BlockSpace):
                obj._linspace: BlockSpace = linspace
            else:
                raise TypeError('linspace is not a BlockSpace object.')
        else:
            start = kwargs['linspace_default_start'] \
                    if 'linspace_default_start' in kwargs.keys() else 0
            stop = kwargs['linspace_default_stop'] \
                    if 'linspace_default_stop' in kwargs.keys() else \
                    start + obj.samples.ncols
            state = kwargs['linspace_default_state'] \
                    if 'linspace_default_state' in kwargs.keys() else \
                    "1"
            obj._linspace: BlockSpace = BlockSpace(start, stop, state)
        return obj

    # Override from_fasta classmethod
    @classmethod
    def from_fasta(cls, path, name=None, comment_parser=None, marker_kw='',
                   **kwargs):
        if not isinstance(marker_kw, str):
            raise TypeError('marker_kw must be a str.')
        return super().from_fasta(
            path, name,comment_parser=comment_parser, keywords=[marker_kw], 
            **kwargs)

    def get_subset(self, samples=None, markers=None, cols=None):
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
        return subset(self, cols, samples=samples, markers=markers)


# TODO: Refactor CatAlignmnet to use mixins
class CatAlignment(Alignment):
    def __init__(self, name, sample_alignment, marker_alignment,
                 linspace=None, subspaces=None, metadata=None, **kwargs):
        super().__init__(name, sample_alignment, marker_alignment, linspace, metadata, **kwargs)
        self._subspaces: OrderedDict = subspaces if subspaces else OrderedDict()

    def subspace(self, name):
        """Returns the linspace of the specified subalignment"""
        return self._subspaces[name]

    @classmethod
    def from_fasta(cls, path, name=None, marker_kw=None, comment_parser=None):
        """Create a CatAlignment object from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file.
        name : str, optional
            Name of the new alignment.
            (default is None, takes the name from the comments
            or uses the filename)
        marker_kw : str, optional
            A sample is considered a marker if this keyword is found
            in the identifier. (default is None, uses the built-in
            comment parser)

        Returns
        -------
        Alignment
            Creates a new CatAlignment object based on the identifiers,
            descriptions, and sequences in the FASTA file.

        """
        if marker_kw is None:
            marker_kw = ''
        # Create alignments
        samples, markers, comment_list = \
            fasta_file_to_basealignments(path, marker_kw)
        comment_parser = \
            comment_parser if comment_parser is not None else parse_cat_comment_list
        kwargs = comment_parser(comment_list)
        if name is not None:
            name = name
        elif 'name' in kwargs.keys():
            name = kwargs['name']
        else:
            name = os.path.basename(path)
        return cls(name, samples, markers, **kwargs)

    def to_fasta(self, path, include_markers=True, 
                 include_headers=True,
                 include_subspaces=True,
                 include_metadata=True):
        """Saves the concatenated alignment as a FASTA-formatted file.
        Some metadata may not be lost.

        Parameters
        ----------
        path : str
            Path to save the alignment to.
        include_markers : bool, optional
            Whether or not to output marker sequences.
            (default is True, include markers in the output)
        include_headers : bool, optional
            Whether or not to output header infomation,
            ie. alignment name and coordinates.
            (default is True, include headers in the output as
            comments)
        include_metadata : bool, optional
            Whether or not to output metadata as comments.
            (default is True, include metadata in the output as
            comments)

        """
        headers_d = {
            'name': str(self.name),
            'coords': '{' + self._linspace.to_block_str() + '}',
        }
        with open(path, 'w') as writer:
            if include_headers:
                for k, v in headers_d.items():
                    print(';{k}\t{v}'.format(k=k, v=v), file=writer)
            if include_subspaces:
                for k, subspace in self._subspaces.items():
                    k = 'subcoords:{}'.format(k)
                    v = '{' + subspace.to_simple_block_str() + '}'
                    print(';{k}\t{v}'.format(k=k, v=v), file=writer)
            if include_metadata:
                for k, v in self.metadata.items():
                    print(';{k}\t{v}'.format(k=k, v=v), file=writer)
            print(self.samples, file=writer)
            if include_markers:
                print(self.markers, file=writer)

    def split_alignment(self):
        """Splits the concatenated alignment into a list of alignments.

        Returns
        -------
        list of Alignment

        """
        aln_list = []
        for name, start, stop in self._linspace.to_list():
            aln = self.get_cols(list(range(start, stop)))
            aln.name = name
            aln._linspace = self._subspaces[name]
            aln.metadata = deepcopy(self.metadata)

            aln_list.append(aln)
        return aln_list

    def __str__(self):
        parts = [
            ';name\t' + str(self.name),
            ';cat_coords\t{' + self._linspace.to_block_str() + '}',
            str(self.samples),
        ]
        if self.markers:
            return '\n'.join(parts + [str(self.markers)])
        return '\n'.join(parts)
