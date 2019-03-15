import pandas as pd

from libalignmentrs.record import Record
from alignmentrs.aln.mixins import serde
from alignmentrs.aln.mixins.tests.mocks import MockData


class MockClass(serde.DictSerdeMixin):
    def __init__(self, records, name=None, index=None, comments=None, row_metadata=None, column_metadata=None, store_history=True, **kwargs):
        self.name = name
        self.index = index
        self.comments = comments
        self.row_metadata = row_metadata
        self.column_metadata = column_metadata
        self.store_history = store_history
        self.kwargs = kwargs
        self.data = MockData()


class TestDictSerdeMixin:

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

    def test_to_dict_with_meta(self):
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test_dict = test_class.to_dict(row_metadata=True, column_metadata=True)
        exp_dict = {
            'name': self.name,
            'data': self.sequences,
            'comments': self.comments,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        }
        assert exp_dict == test_dict, \
            "expected and test dictionaries are not the same: {} != {}".format(
                exp_dict, test_dict
            )

    def test_to_dict_row_meta(self):
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test_dict = test_class.to_dict(row_metadata=True, column_metadata=False)
        exp_dict = {
            'name': self.name,
            'data': self.sequences,
            'comments': self.comments,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': None,
            'column_metadata_index': self.column_metadata.index.to_list(),
        }
        assert exp_dict == test_dict, \
            "expected and test dictionaries are not the same: {} != {}".format(
                exp_dict, test_dict
            )
    
    def test_to_dict_col_meta(self):
        test_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        test_dict = test_class.to_dict(row_metadata=False, column_metadata=True)
        exp_dict = {
            'name': self.name,
            'data': self.sequences,
            'comments': self.comments,
            'row_metadata': None,
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        }
        assert exp_dict == test_dict, \
            "expected and test dictionaries are not the same: {} != {}".format(
                exp_dict, test_dict
            )

    def test_from_dict(self):
        test_dict = {
            'name': self.name,
            'data': self.sequences,
            'comments': self.comments,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        }
        test_class = MockClass.from_dict(test_dict)
        exp_class = MockClass(
            self.records, name=self.name,
            index=self.index,
            comments=self.comments,
            row_metadata=self.row_metadata,
            column_metadata=self.column_metadata,
            store_history=True,
        )
        assert type(exp_class) == type(test_class), \
            "expected and test classes are not the same: {} != {}".format(
                exp_class.__class__.__name__,
                test_class.__class__.__name__,
            )
        # NOTE: Dictionarie values are not compared because special comparison is needed to compare pandas DataFrames 
        assert exp_class.__dict__.keys() == test_class.__dict__.keys(), \
            "expected and test class dictionaries are not the same: {} != {}".format(
                exp_class.__dict__,
                test_class.__dict__,
            )
        assert exp_class.__dir__() == test_class.__dir__(), \
            "expected and test class dir are not the same: {} != {}".format(
                exp_class.__dir__(),
                test_class.__dir__(),
            ) 
