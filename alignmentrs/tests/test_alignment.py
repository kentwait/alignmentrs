import os

import pandas as pd

from alignmentrs.aln import Alignment
from libalignmentrs.alignment import SeqMatrix

def type_error(expected, actual):
    return 'Expected type {}, instead got {}'.format(expected, actual)

def value_error(expected, actual):
    return 'Expected value {}, instead got {}'.format(expected, actual)


class TestMakeData:

    def setup(self):
        # Just initialize the object, dont bother with its data
        self.aln = Alignment(SeqMatrix([]))

    def teardown(self):
        pass

    def test_none(self):
        exp = SeqMatrix([])
        test = self.aln._make_data()
        # Fails on exp == test cmp
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_seqmatrix(self):
        data = SeqMatrix([
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ])
        exp = data
        test = self.aln._make_data(data)
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )
    
    def test_list_str(self):
        exp = SeqMatrix([
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ])
        data = [
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ]
        test = self.aln._make_data(data)
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_tuple_str(self):
        exp = SeqMatrix([
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ])
        data = (
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        )
        test = self.aln._make_data(data)
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_list_list_str(self):
        exp = SeqMatrix([
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ])
        data = [
            list("ATCG"),
            list("ATGG"),
            list("ATCC"),
            list("TAGC"),
        ]
        test = self.aln._make_data(data)
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_list_mixed_str(self):
        exp = SeqMatrix([
            "ATCG",
            "ATGG",
            "ATCC",
            "TAGC",
        ])
        data = [
            list("ATCG"),
            tuple("ATGG"),
            list("ATCC"),
            list("TAGC"),
        ]
        test = self.aln._make_data(data)
        assert exp.data == test.data, \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )


class TestMakeRowMeta:
    def setup(self):
        # Just initialize the object, dont bother with its data
        self.aln = Alignment(SeqMatrix([
            'ATCG',
            'ATGG',
            'ATCC',
            'TAGC',
        ]))

    def teardown(self):
        pass

    def test_none(self):
        exp = pd.DataFrame(None, index=range(4))
        test = self.aln._make_row_meta()
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dataframe(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }, index=range(4))
        data = exp
        test = self.aln._make_row_meta(data)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dict(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        })
        data = {
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }
        test = self.aln._make_row_meta(data)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dict_with_ids(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }, index=[1,2,3,4])
        data = {
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }
        ids = [1,2,3,4]
        test = self.aln._make_row_meta(data, ids=ids)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_with_descriptions_no_ids(self):
        exp = pd.DataFrame({
            'description': [0,1,2,3],
        }, index=[1,2,3,4])
        descriptions = [0,1,2,3]
        test = self.aln._make_row_meta(descriptions=descriptions)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_with_ids_no_descriptions(self):
        exp = pd.DataFrame(None, index=[1,2,3,4])
        ids = [1,2,3,4]
        test = self.aln._make_row_meta(ids=ids)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )


class TestMakeColMeta:
    def setup(self):
        # Just initialize the object, dont bother with its data
        self.aln = Alignment(SeqMatrix([
            'ATCG',
            'ATGG',
            'ATCC',
            'TAGC',
        ]))

    def teardown(self):
        pass

    def test_none(self):
        exp = pd.DataFrame(None, index=range(4))
        test = self.aln._make_col_meta()
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dataframe(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }, index=range(4))
        data = exp
        test = self.aln._make_col_meta(data)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dict(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        })
        data = {
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }
        test = self.aln._make_col_meta(data)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_dict_with_ids(self):
        exp = pd.DataFrame({
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }, index=[1,2,3,4])
        data = {
            'a': [0,1,2,3],
            'b': [0.0,0.1,0.2,0.3],
        }
        ids = [1,2,3,4]
        test = self.aln._make_col_meta(data, ids=ids)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_with_descriptions_no_ids(self):
        exp = pd.DataFrame({
            'description': [0,1,2,3],
        }, index=[1,2,3,4])
        descriptions = [0,1,2,3]
        test = self.aln._make_col_meta(descriptions=descriptions)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )

    def test_with_ids_no_descriptions(self):
        exp = pd.DataFrame(None, index=[1,2,3,4])
        ids = [1,2,3,4]
        test = self.aln._make_col_meta(ids=ids)
        assert all(exp.reset_index() == test.reset_index()), \
            "expected and test SeqMatrix are not the same: {} != {}".format(
                exp, test
            )
