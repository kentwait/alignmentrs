import json
import tempfile

import pandas as pd

from libalignmentrs.record import Record
from alignmentrs.aln.mixins import serde
from alignmentrs.aln.mixins.tests.mocks import MockData


class MockClass(serde.FastaSerdeMixin):
    def __init__(self, records, name=None, index=None, comments=None, row_metadata=None, column_metadata=None, store_history=True, **kwargs):
        self.name = name
        self.index = index
        self.comments = comments
        self.row_metadata = row_metadata
        self.column_metadata = column_metadata
        self.store_history = store_history
        self.kwargs = kwargs
        self.data = MockData()


class TestFastaSerdeMixin:

    def setup(self):
        self.records = [
            Record('test1', 'description1', 'ATGCAT'),
            Record('test2', 'description2', 'ATGGGT'),
            Record('test3', 'description3', 'ATGAAT'),
        ]
        self.name = 'mock_aln'
        self.index = [0,1,2,3,4,5]

        self.comments = {'test_comment': 'testing'}
        self.comments_empty = {}
        self.comments_none = None

        self.row_metadata = pd.DataFrame(
            {
                'description': ['description1', 'description2', 'description3'],
            },
            index=['test1', 'test2', 'test3']
        )
        self.row_metadata_empty = pd.DataFrame({})
        self.row_metadata_none = None

        self.column_metadata = pd.DataFrame(
            {
                'a': [0,1,2,3,4,5],
                'b': [10,11,12,13,14,15],
            },
            index=self.index
        )
        self.column_metdata_empty = pd.DataFrame({})
        self.column_metdata_none = None
    
        self.sequences = [
            'ATGCAT',
            'ATGGGT',
            'ATGAAT',
        ]

        # self.kwargs = kwargs

    def teardown(self):
        pass

    def test_to_fasta(self):
        pass

    def test_from_fasta(self):
        pass

    def test_encode_column_metadata_as_string_str(self):
        exp = "meta|a=[0,1,2,3,4,5]"
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test = test_class._encode_column_metadata_as_string(
            'a', [0,1,2,3,4,5], str
        )
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_encode_column_metadata_as_string_func(self):
        exp = "meta|a=0|1|2|3|4|5"
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        func = lambda lst: '|'.join(list(map(str, lst)))
        test = test_class._encode_column_metadata_as_string(
            'a', [0,1,2,3,4,5], func
        )
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_column_as_string_str(self):
        exp = "meta|a=[0,1,2,3,4,5] meta|b=[10,11,12,13,14,15]"
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test = test_class._column_as_string(['a', 'b'])
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )
    
    def test_column_as_string_func(self):
        exp = "meta|a=[0,1,2,3,4,5] meta|b=10|11|12|13|14|15"
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test = test_class._column_as_string(
            ['a', 'b'],
            encoders={
                'a': str,
                'b': lambda lst: '|'.join(list(map(str, lst))),
            }
        )
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

