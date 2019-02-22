from collections import OrderedDict
import os
import json

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import BaseRecord
from libalignmentrs.position import BlockSpace
from libalignmentrs.position import simple_block_str_to_linspace
from libalignmentrs.readers import fasta_file_to_basealignments


__all__ = ['FastaSerdeMixin', 'DictSerdeMixin', 'JsonSerdeMixin']


class FastaSerdeMixin:
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


class DictSerdeMixin:
    @classmethod
    def from_dict(cls, d):
        records = [BaseRecord(r['id'], r['description'], r['sequence'], 
                              r['chunk_size'])
                   for r in d['alignment']]
        column_metadata = \
            d['column_metadata'] if 'column_metadata' in d.keys() else None
        return cls(d['name'], records, chunk_size=d['chunk_size'],
                   index=d['index'], metadata=d['metadata'],
                   column_metadata=column_metadata)

    def to_dict(self, column_metadata=True):
        d = {
            'name': self.name,
            'index': self._index.to_list(),
            'alignment': [
                {'id': r.id, 'description': r.description,
                 'sequence': r.sequence, 'chunk_size': r.chunk_size}
                for r in self._alignments.records
            ],
            'metadata': self.metadata.to_dict(),
        }
        if column_metadata:
            d['column_metadata'] = self.column_metadata.to_dict()
        return d


class JsonSerdeMixin(DictSerdeMixin):
    @classmethod
    def from_json(cls, path):
        with open(path, 'r') as reader:
            d = json.load(reader)
        return cls.from_dict(d)

    def to_str_json(self, column_metadata=True):
        d = self.to_dict(column_metadata)
        return json.dumps(d)

    def to_json(self, path, column_metadata=True):
        with open(path, 'w') as writer:
            print(self.to_str_json(column_metadata), file=writer)


class NexusSerdeMixin:
    pass


class PhylipSerdeMixin:
    pass


class CsvSerdeMixin:
    @classmethod
    def from_csv(cls, path, delimiter='\t', metadata=True, column_metadata=True):
        # self, name, records, chunk_size: int=1,
        # index=None, metadata: dict=None, column_metadata=None,
        data_path = os.path.abspath(path)

        metadata_d = dict()
        name = None
        chunk_size = 1
        index = None
        records = []
        in_sequence = False
        with open(data_path, 'r') as reader:
            for line in reader.readlines():
                line = line.rstrip()
                if line.startswith('#'):
                    line = line.lstrip()
                    key, value = [string.strip() for string in line.lstrip().split('=')]
                    if key == 'name':
                        name = value
                    elif key == 'chunk_size':
                        chunk_size = value
                    elif key == 'index':
                        # TODO: Add re test for security
                        index = pandas.Index(eval(value))
                    else:
                        metadata[key] = value
                elif line.startswith('id'):
                    in_sequence = True
                    continue
                elif in_sequence is True:
                    sid, desc, seq = line.split(delimiter)
                    records.append(BaseRecord(sid, desc, seq, chunk_size))
        if column_metadata:
            l, r = path.rsplit('.', 1)
            meta_path = l + '.cols.' + r
            column_metadata_d = pandas.read_csv(meta_path, index_col=0, comment='#')
        return cls(name, records, chunk_size=chunk_size,
                   index=index, metadata=metadata_d,
                   column_metadata=column_metadata_d)


    def to_csv(self, path, delimiter='\t', metadata=True, column_metadata=True):
        lines = []
        lines.append('# name = {}'.format(self.name))
        lines.append('# index = {}'.format(self.index.to_list()))
        lines.append('# chunk_size = {}'.format(self._alignment.chunk_size))
        if metadata is True:
            lines += [
                '# {} = {}'.format(k, v) for k, v in metadata.items()
            ]
        lines.append(delimiter.join(['id', 'description', 'sequence']))
        lines += [
            delimiter.join([r.id, r.description.replace('\t', ' '),
                            r.sequence])
            for r in self._alignment.records
        ]
        dirpath = os.path.dirname(path)
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        if not path.endswith('.csv'):
            data_path = path + '.csv'
            meta_path = path + '.cols.csv'
        else:
            data_path = path
            l, r = path.rsplit('.', 1)
            meta_path = l + '.cols.' + r
        with open(data_path, 'w') as writer:
            for line in lines:
                print(line, file=writer)
        if column_metadata is True:
            self._column_metadata.to_csv(meta_path)
