from collections import OrderedDict, ChainMap
from copy import deepcopy
import os
import json
import pickle
import re
import io

import pandas

from libalignmentrs.alignment import SeqMatrix
from libalignmentrs.readers import fasta_to_dict
from alignmentrs.utils import to_intlist


__all__ = [
    'FastaSerdeMixin', 'DictSerdeMixin', 'JsonSerdeMixin', 
    'PickleSerdeMixin', 'CsvSerdeMixin', 'RecordsSerdeMixin',
]


_whitespace_regexp = re.compile(r'\s+')
_column_metadata_string_regexp = re.compile(r'meta\|(\S+)\=(\S+)')

class RecordsSerdeMixin:
    @classmethod
    def from_records(cls, records, name=None, index=None, comments=None, 
                     row_metadata=None, column_metadata=None,
                     store_history=True, **kwargs):
        """Create a new alignment from a list of records.

        Parameters
        ----------
        records : list of Record
        name : str, optional
        index : pandas.Index, optional
            Alingment column position index
        comments : dict, optional
            Comments about the alignment.
        row_metadata : pandas.DataFrame, optional
            Metadata about alignment records.
        column_metadata : pandas.DataFrame, optional
            Metadata for columns in the alignment.
        store_history : bool, optional
            Whether or not to record the actions performed on the alignment.

        Returns
        -------
        Alignment

        """

        return cls(records, name=name, index=index, comments=comments, 
            row_metadata=row_metadata, column_metadata=column_metadata,
            store_history=store_history, **kwargs)

    def to_records(self):
        return [Record(vals[0], vals[1]['description'], self.data.get_row(i))
                for i, vals in enumerate(self.row_metadata.iterrows())]

class FastaSerdeMixin:
    """Adds ability to read/write an Alignment object
    from a FASTA formatted file.
    """
    @classmethod
    def from_fasta(cls, path, name=None, parse_row_metadata=True,  parse_column_metadata=True, store_history=True,column_metadata_decoders=None, **kwargs):
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
        matrix, metadata = fasta_to_dict(path)
        row_meta, col_meta = None, None
        if parse_column_metadata:
            # Parses metadata['descriptions'] and removes parsed info
            pass
        if name is None:
            name = os.path.basename(path)
        return cls(matrix, name,
                   row_metadata=row_meta,
                   # row_ids and row_descriptions are ignored
                   # if row_meta is not None
                   row_ids=metadata['ids'],
                   row_descriptions=metadata['descriptions'],
                   # col_meta is None unless parse_column_metadata is True
                   col_metadata=col_meta,
                   aln_metadata=metadata['comments'],
                   store_history=store_history, **kwargs)

    def to_fasta(self, path=None, include_column_metadata=None, column_metadata_encoders=None, **kwargs):
        """Saves the alignment as a FASTA-formatted file.
        Some metadata may not be lost.

        Parameters
        ----------
        path : str
            Path to save the alignment to.
        include_column_metadata : list of str, optional
            List of keys of columns in column metadata to include
            (default is None, information are not written as FASTA comments
            to ensure maximum compatibility)
        column_metadata_encoders : dict of callable, optional
            Dictionary of functions used to transform values of included
            column metadata.
            Keys are expected to match  specified column names.
            (default is None, all included columns will be transformed using the
            `str` string constructor)

        """
        # Default values if not specified
        if include_column_metadata is None:
            include_column_metadata = []
        if column_metadata_encoders is None:
            column_metadata_encoders = {}

        # Transform col metadata DataFrame into a stringed representation
        # of columns and values.
        col_meta_str = col_metadata_to_str(
            self.column_metadata, include_column_metadata,
            column_metadata_encoders
        )
        # Creates a generator that writes each entry as a string
        # in the FASTA format:
        # >{sid} {desc}
        # {seq}
        info_generator = (
            (vals[0], vals[1]['description'], self.data.data[i]) 
            for i, vals in enumerate(self.row_metadata.iterrows())
        )
        fasta_str = '\n'.join([
            self._entry_formatter(sid, desc, col_meta_str, seq)
            for sid, desc, seq in info_generator
        ])

        # Write the FASTA string to file
        if path is None:
            return fasta_str
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
            print(fasta_str, file=writer)

    @staticmethod
    def _fasta_entry_formatter(sid, desc, col_meta, seq):
        if len(desc) > 0:
            if len(col_meta) > 0:
                return '>{} {} {}\n{}'.format(sid, desc, col_meta, seq)
            return '>{} {}\n{}'.format(sid, desc, seq)

        if len(col_meta) > 0:
            return '>{} {}\n{}'.format(sid, col_meta, seq)
        
        return '>{}\n{}'.format(sid, seq)


class DictSerdeMixin:
    """Adds ability to read/write an Alignment object from a dictionary.
    """
    @classmethod
    def from_dict(cls, d, store_history=True, **kwargs):
        return cls(d['data'],
                   name=d['name'],
                   row_ids=d['row_metadata_index'],
                   row_descriptions=d['row_metadata'],
                   col_ids=d['column_metadata_index'],
                   col_descriptions=d['column_metadata'],
                   aln_metadata=d['alignment_metadata'],
                   store_history=store_history,
                   **kwargs)

    def to_dict(self, row_metadata=True, column_metadata=True):
        d = {
            'name': self.name,
            'data': self.data.data,
            'alignment_metadata': self.alignment_metadata,
        }
        if row_metadata:
            d['row_metadata'] = self.row_metadata.to_dict(orient='list')
            d['row_metadata_index'] = self.row_metadata.index.to_list()
        if column_metadata:
            d['column_metadata'] = self.column_metadata.to_dict(orient='list')
            d['column_metadata_index'] = self.column_metadata.index.to_list()
        return d


class JsonSerdeMixin(DictSerdeMixin):
    @classmethod
    def from_json(cls, path, store_history=True, **kwargs):
        if isinstance(path, io.IOBase):
            d = json.load(path)
        elif isinstance(path, str):
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
        if isinstance(path, io.IOBase):
            obj = pickle.load(path)
        elif isinstance(path, str):
            with open(path, 'rb') as reader:
                obj = pickle.load(reader)
        return obj

    def to_pickle(self, path=None, **kwargs):
        pickled = pickle.dumps(self)
        if path is None:
            return pickled
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'wb') as writer:
            print(pickled, file=writer)

    def __getstate__(self):
        d = {k: v for k, v in self.__dict__.items() if k != 'data'}
        d['data'] = self.data.sequences
        return d

    def __setstate__(self, d):
        d['data'] = from_list(d['data'])
        self.__dict__ = d


class NexusSerdeMixin:
    pass


class PhylipSerdeMixin:
    pass


def col_metadata_to_str(column_metadata, included_keys, encoders=None):
    """Transforms the column metadata DataFrame into a string representation.
    
    Parameters
    ----------
    column_metadata : pandas.DataFrame
        Column metadata
    included_keys : list of str
        List of column metadata column names that to be included in the
        string representation.
    encoders : dict of callable, optional
        Dictionary of functions used to transform column metadata values.
        Keys are expected to match the column names of the column metadata
        DataFrame. (default is None, all columns will be transformed using the
        `str` string constructor)
    
    Returns
    -------
    str
        Column metadata categories and values represented as a string.

    """
    # Creates a tuple generator giving the filtered set of column metadata
    # Each tuple generated consists of the column name and the list of values
    # for that column. 
    included_values = (
        (k, v) for k, v in column_metadata.to_dict(orient='list').items()
        if k in included_keys
    )
    if encoders is None:
        encoders = dict()
    # Creates a list of stringed column metadata where each string is the
    # contains the stringed data of a column metadata category (column)
    # The metadata column is transformed into a string by consuming
    # the `included_values` generator and calling `col_metadata_str_formatter`
    # for each item yielded.
    str_list = [
        col_metadata_str_formatter(
            k, v, encoders[k] if k in encoders.keys() else str)
        for k, v in included_values
    ]
    # Each column's string representation is separated by a whitespace
    return ' '.join(str_list)

def col_metadata_str_formatter(key, value, encoder: callable):
    """Returns the string representation of a column metadata category.
    
    Parameters
    ----------
    key : str
        Name of column metadata category (column name in the DataFrame).
    value : list
        Column metadata category values.
    encoder : callable
        Function used to transform the list of values into a string.
    
    Returns
    -------
    str
        String representation of the column metadata category and its values.

    """
    return 'c|{}={}'.format(key, _whitespace_regexp.sub('', encoder(value)))

def make_col_meta_dict(description, decoders):
    matches = _column_metadata_string_regexp.findall(description)
    return {
        k: (decoders[k](v) if k in decoders.keys() else eval(v))
        for k, v in matches
    }
