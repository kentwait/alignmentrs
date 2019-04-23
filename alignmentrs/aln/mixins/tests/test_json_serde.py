import json
import tempfile

import pandas as pd

from libalignmentrs.record import Record
from alignmentrs.aln.mixins import serde
from alignmentrs.aln.mixins.tests.mocks import MockData


class MockAlignment(serde.JsonSerdeMixin):
    def __init__(self, matrix, name=None,
                 row_metadata=None, col_metadata=None,
                 row_ids:list=None, row_descriptions:list=None,
                 col_ids:list=None, col_descriptions:list=None,
                 aln_metadata:dict=None, store_history=True,
                 **kwargs):
        self.name = name
        self.alignment_metadata = aln_metadata
        if row_metadata is None:
            self.row_metadata = pd.DataFrame(row_descriptions, index=row_ids)
        else:
            self.row_metadata = row_metadata
        if col_metadata is None:
            self.column_metadata = pd.DataFrame(col_descriptions, index=col_ids)
        else:
            self.column_metadata = col_metadata
        self.store_history = store_history
        self.kwargs = kwargs
        self.data = MockData()

class TestJsonSerdeMixin:

    def setup(self):
        self.name = 'mock_aln'
        self.alignment_metadata = {'comment1': 'This is a comment.'}
        self.row_metadata = pd.DataFrame(
            {'description': ['','d2','desc3']}, index=['test1','test2','test3'])
        self.column_metadata = pd.DataFrame({'a': [1,2,3,4,5,6]}, index=range(6))
        self.store_history = True
        # self.kwargs = kwargs
        self.matrix = MockData()
        self.test_aln = MockAlignment(
            self.matrix, name=self.name,
            row_metadata=self.row_metadata,
            col_metadata=self.column_metadata,
            aln_metadata=self.alignment_metadata,
            store_history=self.store_history,
        )

    def teardown(self):
        pass

    def test_to_json_with_row_col_meta(self):
        test_json = self.test_aln.to_json(
            row_metadata=True, column_metadata=True)
        exp_json = json.dumps({
            'name': self.name,
            'data': self.matrix.data,
            'alignment_metadata': self.alignment_metadata,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        })
        assert exp_json == test_json, \
            "expected and test json are not the same: {} != {}".format(
                exp_json, test_json
            )

    def test_to_dict_with_row_meta(self):
        test_json = self.test_aln.to_json(
            row_metadata=True, column_metadata=False)
        exp_json = json.dumps({
            'name': self.name,
            'data': self.matrix.data,
            'alignment_metadata': self.alignment_metadata,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
        })
        assert exp_json == test_json, \
            "expected and test dictionaries are not the same: {} != {}".format(
                exp_json, test_json
            )
    
    def test_to_dict_with_col_meta(self):
        test_json = self.test_aln.to_json(
            row_metadata=False, column_metadata=True)
        exp_json = json.dumps({
            'name': self.name,
            'data': self.matrix.data,
            'alignment_metadata': self.alignment_metadata,
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        })
        assert exp_json == test_json, \
            "expected and test dictionaries are not the same: {} != {}".format(
                exp_json, test_json
            )

    def test_from_json(self):
        exp_json = json.dumps({
            'name': self.name,
            'data': self.matrix.data,
            'alignment_metadata': self.alignment_metadata,
            'row_metadata': self.row_metadata.to_dict(orient='list'),
            'row_metadata_index': self.row_metadata.index.to_list(),
            'column_metadata': self.column_metadata.to_dict(orient='list'),
            'column_metadata_index': self.column_metadata.index.to_list(),
        })
        with tempfile.TemporaryFile(mode='r+') as f:
            f.write(exp_json)
            f.seek(0)
            test_class = MockAlignment.from_json(f)

        exp_class = self.test_aln
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
