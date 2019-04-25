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
    'col_metadata_to_str', 'col_metadata_str_formatter',
]


_whitespace_regexp = re.compile(r'\s+')
_column_metadata_string_regexp = re.compile(r'meta\|(\S+)\=(\S+)')


class FastaSerdeMixin:
    """Adds ability to read/write an Alignment object
    from a FASTA formatted file.
    """
    @classmethod
    def from_fasta(cls, path, name=None, parse_row_metadata=True,  parse_description=True, column_metadata_decoders=None, column_metadata_regexp='c\|([A-Za-z0-9\s\.]+)=(\[[A-Za-z0-9\.\s,\"\']+\])', column_index_regexp='ci\|([A-Za-z0-9\s\.]+)=(\[[A-Za-z0-9\.\s,\"\']+\])', store_history=True, **kwargs):
        """Create an Alignment object from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file.
        name : str, optional
            Name of the new alignment.
            (default is None, takes the name from the comments
            or uses the filename)
        parse_description : function, optional
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
        if parse_description:
            # Parses metadata['descriptions'] and removes parsed info
            offset = 0
            match_locations = []
            col_d = {}
            col_idx = None
            # Parses column index
            match = re.search(column_index_regexp, metadata['descriptions'][0])
            if match:
                key, value = match.groups()
                try:
                    value = eval(value)
                except SyntaxError:
                    raise ValueError('Cannot construct Alignment from the given FASTA file: column index is malformed'.format(key))
                # Put key-value pair into the dictionary
                col_idx = value

            # Parses column metadata
            for match in re.finditer(column_metadata_regexp,
                                     metadata['descriptions'][0]):
                match_locations.append(match.span())
                key, value = match.groups()
                # Convert text into a list using eval
                # This is DANGEROUS and could open to exploitation.
                # TODO: Add a prelimenary regex check to lessen vulnerability
                try:
                    value = eval(value)
                except SyntaxError:
                    raise ValueError('Cannot construct Alignment from the given FASTA file: column metadata {} is malformed'.format(key))
                # Put key-value pair into the dictionary
                col_d[key] = value
            # Constructs column metadata DataFrame from dictionary and index
            if (col_idx is not None) or col_d:
                col_meta = pandas.DataFrame(col_d, index=col_idx)

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

    def to_fasta(self, path=None, include_column_metadata=None, column_metadata_encoders=None, column_metadata_template='c|{}={}', **kwargs):
        """Saves the alignment as a FASTA-formatted file.
        Some metadata may not be lost.

        Parameters
        ----------
        path : str, optional
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
            column_metadata_encoders, column_metadata_template
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
            self._fasta_entry_formatter(*params, col_meta_str)
                if i == 0 else
                self._fasta_entry_formatter(*params, '')
            for i, params in enumerate(info_generator)
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
    def _fasta_entry_formatter(sid, desc, seq, col_meta):
        # Formats the ID, description, stringed metadata, and sequence
        # to follow the FASTA format.
        # There are 4 possible scenarios, note that the identifier string
        # is always expected to have a non-empty value:
        # - all information exists
        # - col_meta is an empty string
        # - description is an empty string
        # - both col_meta and description are empty strings

        # Checks if ID is empty
        if len(sid) < 1:
            raise ValueError('Cannot create FASTA file: identifier string cannot be empty.')
        # Description is not empty
        if len(desc) > 0:
            if len(col_meta) > 0:
                return '>{} {} {}\n{}'.format(sid, desc, col_meta, seq)
            return '>{} {}\n{}'.format(sid, desc, seq)
        # Description is empty but column metadata is not empty
        if len(col_meta) > 0:
            return '>{} {}\n{}'.format(sid, col_meta, seq)
        # Decription and column metadata are empty
        return '>{}\n{}'.format(sid, seq)


class DictSerdeMixin:
    """Adds ability to read/write an Alignment object from a dictionary.
    """
    @classmethod
    def from_dict(cls, d, store_history=True, **kwargs):
        """Creates an Alignment object from a dictionary.

        Parameters
        ----------
        d : dict
            Dictionary containing the alignment information and relevant
            metadata.

        Returns
        -------
        Alignment

        """
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
        """Returns the dictionary representation of the alignment.
        Contents of the dictionary use builtin types to maximize
        compatibility.

        Parameters
        ----------
        row_metadata : bool, optional
            Whether or not to include row metadata information. (default is True, row metadata is included)
        column_metadata : bool, optional
            Whether or not to include column metadata information. (default is True, column metadata is included)

        Returns
        -------
        dict

        """
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
    """Adds ability to read/write an Alignment object from a JSON file.

    The underlying organization of the JSON encoding is based on the dictionary 
    created using the DictSerdeMixin dictionary mixin.
    """
    @classmethod
    def from_json(cls, path, store_history=True, **kwargs):
        """Create an alignment from a JSON file.
        
        Parameters
        ----------
        path : io.IOBase or str
            File stream using a file handler or a string to the path.
        
        Returns
        -------
        Alignment

        """
        if isinstance(path, io.IOBase):
            # json.load requires a file handler to read the file.
            # io.IOBase is the abstract base class for all kinds of
            # I/O file streaming.
            d = json.load(path)
        elif isinstance(path, str):
            # The user can also input the path where the json file is located.
            # To handle this, the path will have to be opened as a file handler
            with open(path, 'r') as reader:
                d = json.load(reader)
        # JSON structure is based on to_dict and so it is dependent
        # on DictSerdeMixin.
        return cls.from_dict(d, store_history=store_history, **kwargs)

    def to_json(self, path=None, row_metadata=True, column_metadata=True):
        """Saves the alignment as a JSON file.

        Parameters
        ----------
        path : str, optional
            Path to save the alignment to.
        row_metadata : bool, optional
            Whether or not to include row metadata information. (default is True, row metadata is included)
        column_metadata : bool, optional
            Whether or not to include column metadata information. (default is True, column metadata is included)

        Returns
        -------
        str
            If path is None, returns the JSON-formatted text as a string.

        """
        # to_json uses to_dict to transform the alignment data
        # into a representation that uses only builtins to maximize
        # compatibility.
        # The resulting dictionary is encoded into JSON.
        d = self.to_dict(
            row_metadata=row_metadata,
            column_metadata=column_metadata)
        json_str = json.dumps(d)
        # If the save path is not specified, the encoded JSON text
        # is returned as a string
        if path is None:
            return json_str
        # If the save path is specified, the JSON encoded text
        # is written as a text file.
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'w') as writer:
            print(json_str, file=writer)


class PickleSerdeMixin:
    """Adds ability to pickle/unpickle an Alignment object.
    """
    @classmethod
    def from_pickle(cls, path, store_history=True, **kwargs):
        """Converts a pickled alignment back into an Alignment object.

        Parameters
        ----------
        path : io.IOBase or str
            File handler for the pickled alignment or a string to the path.

        Returns
        -------
        Alignment

        """
        if isinstance(path, io.IOBase):
            obj = pickle.load(path)
        elif isinstance(path, str):
            with open(path, 'rb') as reader:
                obj = pickle.load(reader)
        return obj

    def to_pickle(self, path=None, **kwargs):
        """Pickles the current alignment.

        Parameters
        ----------
        path : str, optional
            Path to save the alignment to.

        Returns
        -------
        bytestring
            If path is None, returns the bytestring representation of the
            pickled alignment.

        """
        # Pickles the alignment. Underlying methods that do the pickling are
        # __getstate__ and __setstate__.
        pickled = pickle.dumps(self)
        # If path is not provided, the bytestring of the pickle is returned.
        if path is None:
            return pickled
        # If path is provided, the pickle is written to file.
        dirpath = os.path.dirname(os.path.abspath(path))
        if not os.path.isdir(dirpath):
            raise OSError('{} does not exist'.format(dirpath))
        with open(path, 'wb') as writer:
            print(pickled, file=writer)

    def __getstate__(self):
        # This method gets called when the Alignment object
        # is being pickled.
        d = {k: v for k, v in self.__dict__.items() if k != 'data'}
        d['data'] = self.data.data
        return d

    def __setstate__(self, d):
        # This method gets called when the pickled object
        # is being unpickled back into an Alignment object.
        d['data'] = SeqMatrix(d['data'])
        self.__dict__ = d


class NexusSerdeMixin:
    pass


class PhylipSerdeMixin:
    pass


def col_metadata_to_str(column_metadata, included_keys, encoders=None, template='c|{}={}', index_template='ci|{}={}'):
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
    template : str, optional
        Template used for formatting the string. Template should have 2
        slots for the key and the column value.
    
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
            k, v, encoders[k] if k in encoders.keys() else None, template)
        for k, v in included_values
    ]
    str_index = [col_metadata_str_formatter(
        'index', column_metadata.index.to_list(),
        encoders['index'] if 'index' in encoders.keys() else None, 
        index_template)
    ]
    # Each column's string representation is separated by a whitespace
    return ' '.join(str_index + str_list)

def col_metadata_str_formatter(key, value, encoder:callable=None, template='c|{}={}'):
    """Returns the string representation of a column metadata category.
    
    Parameters
    ----------
    key : str
        Name of column metadata category (column name in the DataFrame).
    value : list
        Column metadata category values.
    encoder : callable
        Function used to transform the list of values into a string.
    template : str, optional
        Template used for formatting the string. Template should have 2
        slots for the key and the column value.

    Returns
    -------
    str
        String representation of the column metadata category and its values.

    """
    if encoder is None:
        encoder = lambda x: _whitespace_regexp.sub('', str(x))
    return template.format(key, encoder(value))

def make_col_meta_dict(description, decoders):
    matches = _column_metadata_string_regexp.findall(description)
    return {
        k: (decoders[k](v) if k in decoders.keys() else eval(v))
        for k, v in matches
    }
