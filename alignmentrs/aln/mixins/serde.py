from collections import OrderedDict
import os
import json
import pickle

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import BaseRecord
from libalignmentrs.readers import fasta_to_records


__all__ = [
    'FastaSerdeMixin', 'DictSerdeMixin', 'JsonSerdeMixin', 
    'PickleSerdeMixin', 'CsvSerdeMixin'
]


class FastaSerdeMixin:
    @classmethod
    def from_fasta(cls, path, name=None, chunk_size=1, column_metadata=None):
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
        records, _ = fasta_to_records(path, chunk_size)
        if column_metadata is not None and isinstance(column_metadata, dict):
            col_meta = []
            keep_records = []
            for i, record in enumerate(records):
                if record.id.startswith('meta|'):
                    _, k = record.id.lsplit('|', 1)
                    if k in column_metadata.keys():
                        func = column_metadata[k]
                        data = list(map(func, record.sequence))
                        col_meta.append({k: data})
                    else:
                        keep_records.append(i)
                else:
                    keep_records.append(i)
            records = records[keep_records]

        return cls(name, records, chunk_size=chunk_size)

    def to_fasta(self, path=None, column_metadata=None):
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
        fasta_str = '\n'.join([str(rec) for rec in self._alignments.records])
        if column_metadata is not None and isinstance(column_metadata, list):
            meta_str = '\n'.join([
                '>meta|{}\n{}'.format(
                    k, 
                    ''.join(list(map(str, self.column_metadata[k].to_list())))
                )
                for k in column_metadata
            ])
            fasta_str += column_metadata
        if path is None:
            return fasta_str
        dirpath = os.path.dirname(path)
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
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

    def to_json(self, path=None, column_metadata=True):
        d = self.to_dict(column_metadata)
        json_str = json.dumps(d)
        if path is None:
            return json_str
        dirpath = os.path.dirname(path)
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
            print(json_str, file=writer)


class PickleSerdeMixin(DictSerdeMixin):
    @classmethod
    def from_pickle(cls, path):
        with open(path, 'rb') as reader:
            d = pickle.load(reader)
        return cls.from_dict(d)

    def to_pickle(self, path, column_metadata=True):
        d = self.to_dict(column_metadata)
        dirpath = os.path.dirname(path)
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'wb') as writer:
            pickle.dump(d, writer)


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


    def to_csv(self, path=None, delimiter='\t', metadata=True, column_metadata=True):
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
        csv_str = '\n'.join(lines)
        if path is None:
            return csv_str
        dirpath = os.path.dirname(path)
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(data_path, 'w') as writer:
            print(csv_str, file=writer)
        if column_metadata is True:
            self._column_metadata.to_csv(meta_path)
