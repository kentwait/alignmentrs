from collections import OrderedDict, ChainMap
import os
import json
import pickle
import re

import pandas

from libalignmentrs.alignment import BaseAlignment
from libalignmentrs.record import Record
from libalignmentrs.readers import fasta_to_records


__all__ = [
    'FastaSerdeMixin', 'DictSerdeMixin', 'JsonSerdeMixin', 
    'PickleSerdeMixin', 'CsvSerdeMixin', 'RecordsSerdeMixin',
]


class RecordsSerdeMixin:
    @classmethod
    def from_records(cls, records, name=None, index=None, comments=None, 
                     row_metadata=None, column_metadata=None,
                     store_history=True, **kwargs):
        cls(records, name=name, index=index, comments=comments, 
            row_metadata=row_metadata, column_metadata=column_metadata,
            store_history=store_history, **kwargs)

    def to_records(self):
        return [Record(vals[0], vals[1]['description'], self.data.get_row(i))
                for i, vals in enumerate(self.row_metadata.iterrows())]

class FastaSerdeMixin:
    @classmethod
    def from_fasta(cls, path, name=None, parse_row_metadata=True,  
                   parse_column_metadata=True, store_history=True, **kwargs):
        """Create an Alignment object from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file.
        name : str, optional
            Name of the new alignment.
            (default is None, takes the name from the comments
            or uses the filename)
        parse_column_metadata : function, optional
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
        records, _ = fasta_to_records(path)
        row_meta, col_meta = [], []
        # row_meta_regexp = re.compile(r'[^(meta\|)](\S+)\=(\S+)')
        col_meta_regexp = re.compile(r'meta\|(\S+)\=(\S+)')
        # if parse_row_metadata:
        #     for record in records:
        #         matches = dict(row_meta_regexp.findall(record.description))
        #         row_meta.append(matches)
        if parse_column_metadata:
            for record in records:
                matches = dict(col_meta_regexp.findall(record.description))
                col_meta.append(matches)
        # row_meta = {k: eval(v) for k, v in dict(ChainMap(*row_meta)).items()}
        col_meta = {k: eval(v) for k, v in dict(ChainMap(*col_meta)).items()}
        if name is None:
            name = os.path.basename(path)
        return cls(records, name=name,
                #    row_metadata=row_meta,
                   column_metadata=col_meta,
                   store_history=store_history, **kwargs)

    def to_fasta(self, path=None, include_column_metadata=None):
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
        include_column_metadata : list, optional
            List of keys of columns in column metadata to include
            (default is None, information are not written as FASTA comments
            to ensure maximum compatibility)

        """
        sp_sub = re.compile(r'\s+')
        id_desc_seq = ((vals[0], vals[1]['description'], self.data.get_row(i)) 
                       for i, vals in enumerate(self.row_metadata.iterrows()))
        col_meta = ' '.join([
            'meta|{}={}'.format(k, sp_sub.sub('', str(v)))
            for k, v in self.column_metadata.to_dict(orient='list').items()
            if k in include_column_metadata
        ])
        fasta_str = '\n'.join([
            self._record_formatter(sid, desc, seq, col_meta)
            for sid, desc, seq in id_desc_seq
        ])
        # TODO: Add ability to add row metedata thats not the description
        if path is None:
            return fasta_str
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
            print(fasta_str, file=writer)

    @staticmethod
    def _record_formatter(sid, desc, seq, col_meta):
        if len(desc) > 0 and len(col_meta) > 0:
            return '>{} {} {}\n{}'.format(sid, desc, col_meta, seq)
        elif len(desc) == 0 and len(col_meta) > 0:
            return '>{} {}\n{}'.format(sid, col_meta, seq)
        elif len(desc) > 0 and len(col_meta) == 0:
            return '>{} {}\n{}'.format(sid, desc, seq)
        return '>{}\n{}'.format(sid, seq)

class DictSerdeMixin:
    @classmethod
    def from_dict(cls, d, store_history=True, **kwargs):
        name = d['name']
        records = [Record(r['id'], r['description'], r['sequence'])
                   for r in d['data']]
        row_metadata = d['row_metadata']
        column_metadata = d['column_metadata']
        index = d['index']
        return cls(records, name=name, index=index,
                   row_metadata=row_metadata,
                   column_metadata=column_metadata,
                   store_history=store_history,
                   **kwargs)

    def to_dict(self, row_metadata=True, column_metadata=True):
        d = {
            'name': self.name,
            'data': self.data.sequences,
            'comments': self.comments,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'index': self.column_metadata.index.to_list(),
        }
        # TODO: Store history
        return d


class JsonSerdeMixin(DictSerdeMixin):
    @classmethod
    def from_json(cls, path, store_history=True, **kwargs):
        with open(path, 'r') as reader:
            d = json.load(reader)
        return cls.from_dict(d, store_history=store_history,
                             **kwargs)

    def to_json(self, path=None, row_metadata=True, column_metadata=True):
        d = self.to_dict(
            row_metadata=row_metadata,
            column_metadata=column_metadata)
        json_str = json.dumps(d)
        if path is None:
            return json_str
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
            print(json_str, file=writer)


class PickleSerdeMixin(DictSerdeMixin):
    @classmethod
    def from_pickle(cls, path, store_history=True, **kwargs):
        with open(path, 'rb') as reader:
            d = pickle.load(reader)
        return cls.from_dict(d, store_history=store_history,
                             **kwargs)

    def to_pickle(self, path, column_metadata=True):
        d = self.to_dict(column_metadata)
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'wb') as writer:
            pickle.dump(d, writer)


class NexusSerdeMixin:
    pass


class PhylipSerdeMixin:
    pass


# class CsvSerdeMixin:
#     @classmethod
#     def from_csv(cls, path, delimiter='\t', metadata=True,
#                  column_metadata=True, store_history=True,
#                  **kwargs):
#         # self, name, records, chunk_size: int=1,
#         # index=None, metadata: dict=None, column_metadata=None,
#         data_path = os.path.abspath(path)

#         metadata_d = dict()
#         name = None
#         chunk_size = 1
#         index = None
#         records = []
#         in_sequence = False
#         with open(data_path, 'r') as reader:
#             for line in reader.readlines():
#                 line = line.rstrip()
#                 if line.startswith('#'):
#                     line = line.lstrip()
#                     key, value = [string.strip() for string in line.lstrip().split('=')]
#                     if key == 'name':
#                         name = value
#                     elif key == 'chunk_size':
#                         chunk_size = value
#                     elif key == 'index':
#                         # TODO: Add re test for security
#                         index = pandas.Index(eval(value))
#                     else:
#                         metadata[key] = value
#                 elif line.startswith('id'):
#                     in_sequence = True
#                     continue
#                 elif in_sequence is True:
#                     sid, desc, seq = line.split(delimiter)
#                     records.append(Record(sid, desc, seq, chunk_size))
#         if column_metadata:
#             l, r = path.rsplit('.', 1)
#             meta_path = l + '.cols.' + r
#             column_metadata_d = pandas.read_csv(meta_path, index_col=0, comment='#')
#         return cls(name, records, chunk_size=chunk_size,
#                    index=index, metadata=metadata_d,
#                    column_metadata=column_metadata_d, store_history=store_history,
#                    **kwargs)


#     def to_csv(self, path=None, delimiter='\t', metadata=True, column_metadata=True):
#         lines = []
#         lines.append('# name = {}'.format(self.name))
#         lines.append('# index = {}'.format(self.index.to_list()))
#         lines.append('# chunk_size = {}'.format(self._alignment.chunk_size))
#         if metadata is True:
#             lines += [
#                 '# {} = {}'.format(k, v) for k, v in metadata.items()
#             ]
#         lines.append(delimiter.join(['id', 'description', 'sequence']))
#         lines += [
#             delimiter.join([r.id, r.description.replace('\t', ' '),
#                             r.sequence])
#             for r in self._alignment.records
#         ]
#         dirpath = os.path.dirname(os.path.abspath(path))
#         if not os.path.isdir(dirpath):
#             raise OSError('{} does not exist'.format(dirpath))
#         if not path.endswith('.csv'):
#             data_path = path + '.csv'
#             meta_path = path + '.cols.csv'
#         else:
#             data_path = path
#             l, r = path.rsplit('.', 1)
#             meta_path = l + '.cols.' + r
#         csv_str = '\n'.join(lines)
#         if path is None:
#             return csv_str
#         dirpath = os.path.dirname(os.path.abspath(path))
#         if not os.path.isdir(dirpath):
#             raise OSError('{} does not exist'.format(dirpath))
#         with open(data_path, 'w') as writer:
#             print(csv_str, file=writer)
#         if column_metadata is True:
#             self.column_metadata.to_csv(meta_path)
