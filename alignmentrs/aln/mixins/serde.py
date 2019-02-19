from collections import OrderedDict
import os

from libalignmentrs.readers import fasta_file_to_basealignments
from libalignmentrs.position import BlockSpace


__all__ = ['FastaSerde']


class FastaSerde:
    @classmethod
    def from_fasta(cls, path, name=None, comment_parser=None, **kwargs):
        """Create an Alignment object from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file.
        name : str, optional
            Name of the new alignment.
            (default is None, takes the name from the comments
            or uses the filename)
        comment_parser : function, optional
            Function that takes a list of comment lines as input
            and outputs a dictionary that organizes comments into
            keys and values. (default is None, lines starting with 
            a semicolon ";" are ignored.)

        Returns
        -------
        Alignment
            Creates a new Alignment object based on the identifiers,
            descriptions, and sequences in the FASTA file.

        """
        if 'keywords' in kwargs.keys():
            balns, comment_list = \
                fasta_file_to_basealignments(path, kwargs['keywords'])
        else:
            balns, comment_list = \
                fasta_file_to_basealignments(path, [])

        # Parse comment if parser is provided
        # Ignore otherwise
        if comment_parser is not None and callable(comment_parser):
            metadata = comment_parser(comment_list)
        else:
            metadata = OrderedDict()

        # Name the alignment object
        if name is not None:
            name = name
        elif 'name' in metadata.keys():
            name = metadata['name']
        else:
            name = os.path.basename(path)

        # Check if linspace in metadata
        if 'linspace' in metadata.keys() and \
            isinstance(metadata['linspace'], BlockSpace):
            linspace = metadata['linspace']
        else:
            linspace = None

        # Removes name and linspace from metadata if present
        if metadata:
            for k in ['name', 'linspace']:
                del metadata[k]
        
        # Create a new Alignment object
        return cls(name, *balns, linspace=linspace, metadata=metadata)


    def to_fasta_str(self, include_info=False, include_metadata=False):
        parts = []
        if include_info:
            parts += [
                ';name\t' + str(self.name),
                ';coords\t{' + self._linspace.to_simple_block_str() + '}',
            ]
        if include_metadata:
            parts += [';{}\t{}'.format(k, v) for k, v in self.metadata.items()]
        for member in self.__class__.members:
            parts.append(str(self.__getattribute__(member)))
        return '\n'.join(parts)

    def to_fasta(self, path, include_info=False, include_metadata=False):
        """Saves the alignment as a FASTA-formatted file.
        Some metadata may not be lost.

        Parameters
        ----------
        path : str
            Path to save the alignment to.
        include_info : bool, optional
            Whether or not to output alignment infomation,
            ie. alignment name and coordinates.
            (default is False, information are not written as FASTA comments
            to ensure maximum compatibility)
        include_metadata : bool, optional
            Whether or not to output metadata as comments.
            (default is False, information are not written as FASTA comments
            to ensure maximum compatibility)

        """
        with open(path, 'w') as writer:
            fasta_str = self.to_fasta_str(
                include_info=include_info, include_metadata=include_metadata)
            print(fasta_str, file=writer)

