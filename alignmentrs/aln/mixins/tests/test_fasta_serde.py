import json
import tempfile

import pandas as pd

from libalignmentrs.record import Record
from alignmentrs.aln.mixins import serde
from alignmentrs.aln.mixins.serde import (
    col_metadata_str_formatter, col_metadata_to_str)
from alignmentrs.aln.mixins.tests.mocks import MockData


class MockAlignment(serde.FastaSerdeMixin):
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

class TestFastaSerdeMixin:

    def setup(self):
        self.name = 'mock_aln'
        self.alignment_metadata = {'comment1': 'This is a comment.'}
        self.row_metadata = pd.DataFrame(
            {'description': ['test description','','desc3']}, index=['test1','test2','test3'])
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

    def test_to_fasta(self):
        exp = '>test1 test description ci|index=[0,1,2,3,4,5] c|a=[1,2,3,4,5,6]\nATGCAT\n>test2\nATGGGT\n>test3 desc3\nATGAAT'
        test = self.test_aln.to_fasta(include_column_metadata=['a'])
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_from_fasta(self):
        exp_fasta = '>test1 test description ci|index=[0,1,2,3,4,5] c|a=[1,2,3,4,5,6]\nATGCAT\n>test2\nATGGGT\n>test3 desc3\nATGAAT'
        with tempfile.NamedTemporaryFile(mode='r+') as f:
            f.write(exp_fasta)
            f.seek(0)
            test_class = MockAlignment.from_fasta(f.name, name='mock_aln')
        
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

    def test_fasta_entry_formatter_with_desc_meta(self):
        exp = '>test1 this is a description meta|a=[0,1,2,3,4,5]\nATGCGC'
        test = self.test_aln._fasta_entry_formatter(
            'test1', 'this is a description', 'ATGCGC', 'meta|a=[0,1,2,3,4,5]')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_fasta_entry_formatter_with_desc(self):
        exp = '>test1 this is a description\nATGCGC'
        test = self.test_aln._fasta_entry_formatter(
            'test1', 'this is a description', 'ATGCGC', '')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_fasta_entry_formatter_with_meta(self):
        exp = '>test1 meta|a=[0,1,2,3,4,5]\nATGCGC'
        test = self.test_aln._fasta_entry_formatter(
            'test1', '', 'ATGCGC', 'meta|a=[0,1,2,3,4,5]')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_fasta_entry_formatter_no_desc_meta(self):
        exp = '>test1\nATGCGC'
        test = self.test_aln._fasta_entry_formatter(
            'test1', '', 'ATGCGC', '')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )


class TestColMetadataToStr:
    def setup(self):
        self.df = pd.DataFrame({
            'data1': list(range(0,5)),
            'data2': list(map(lambda x: x/100, range(0,5))),
            'd3': list(map(lambda x: x**2, range(0,5))),
        })

    def teardown(self):
        pass

    def test_all_keys_no_encoder(self):
        exp = 'ci|index=[0,1,2,3,4] c|data1=[0,1,2,3,4] c|data2=[0.0,0.01,0.02,0.03,0.04] c|d3=[0,1,4,9,16]'
        test = col_metadata_to_str(self.df, ['data1', 'data2', 'd3'])
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )
    
    def test_all_keys_no_encoder_custom_template(self):
        exp = 'ci|index=[0,1,2,3,4] (col|data1=[0,1,2,3,4]) (col|data2=[0.0,0.01,0.02,0.03,0.04]) (col|d3=[0,1,4,9,16])'
        test = col_metadata_to_str(self.df, ['data1', 'data2', 'd3'], template='(col|{}={})')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )


    def test_all_keys_no_encoder_custom_index_template(self):
        exp = 'COLIDX|index=[0,1,2,3,4] c|data1=[0,1,2,3,4] c|data2=[0.0,0.01,0.02,0.03,0.04] c|d3=[0,1,4,9,16]'
        test = col_metadata_to_str(self.df, ['data1', 'data2', 'd3'], index_template='COLIDX|{}={}')
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_all_keys_custom_encoder(self):
        exp = 'ci|index=[0,1,2,3,4] c|data1=[0, 1, 2, 3, 4] c|data2=[0.0,0.01,0.02,0.03,0.04] c|d3=[0,1,4,9,16]'
        encoders = {
            'data1': lambda x: str(x),
            'data2': lambda x: str(x).replace(' ', ''),
            'd3': None,
        }
        test = col_metadata_to_str(
            self.df, ['data1', 'data2', 'd3'], encoders=encoders)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_some_keys_custom_encoder(self):
        exp = 'ci|index=[0,1,2,3,4] c|data1=[0, 1, 2, 3, 4] c|d3=[0,1,4,9,16]'
        encoders = {
            'data1': lambda x: str(x),
            # extra encoder become data2 column is not included
            'data2': lambda x: str(x).replace(' ', ''),  
            'd3': None,
        }
        test = col_metadata_to_str(
            self.df, ['data1', 'd3'], encoders=encoders)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_some_keys_custom_index_encoder(self):
        exp = 'ci|index=[0, 1, 2, 3, 4] c|data1=[0,1,2,3,4] c|d3=[0,1,4,9,16]'
        encoders = {
            'index': str,
            'data1': None,
            'd3': None,
        }
        test = col_metadata_to_str(
            self.df, ['data1', 'd3'], encoders=encoders)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )


class TestColMetadataStrFormatter:
    def setup(self):
        self.df = pd.DataFrame({
            'data1': list(range(0,5)),
            'data2': list(map(lambda x: x/100, range(0,5))),
            'd3': list(map(lambda x: x**2, range(0,5))),
        })
        self.key = 'data1'
        self.value = list(range(0,5))

    def teardown(self):
        pass

    def test_defaults(self):
        exp = 'c|data1=[0,1,2,3,4]'
        test = col_metadata_str_formatter(self.key, self.value)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )
    
    def test_custom_encoder(self):
        exp = 'c|data1=[0, 1, 2, 3, 4]'
        encoder = str
        test = col_metadata_str_formatter(self.key, self.value, encoder=encoder)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )

    def test_custom_template(self):
        exp = '(col|data1=[0,1,2,3,4])'
        template = '(col|{}={})'
        test = col_metadata_str_formatter(self.key, self.value, template=template)
        assert exp == test, \
            "expected and test strings are not the same: {} != {}".format(
                exp, test
            )
