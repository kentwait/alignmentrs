import os

import pandas as pd

from alignmentrs.aln import Alignment
from libalignmentrs.alignment import SeqMatrix

def type_error(expected, actual):
    return 'Expected type {}, instead got {}'.format(expected, actual)

def value_error(expected, actual):
    return 'Expected value {}, instead got {}'.format(expected, actual)


class TestMakeRowMeta():
    def setUp(self):
        self.matrix = SeqMatrix([
            "atcg",
            "atgg",
            "atcc",
            "tagc",
        ])
        self.aln = Alignment(self.matrix)

    def tearDown(self):
        pass

    def test_defaults(self):
        exp = pd.DataFrame([], index=range(4))
        res = self.aln._make_row_meta()
        assert exp == res
