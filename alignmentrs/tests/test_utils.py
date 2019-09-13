""" Unit test for utils functions. """

import tempfile
from nose.tools import *

from alignmentrs.utils import fasta_file_to_lists, alignment_file_to_lists

class TestFastaFileReader:
    """ Unit tests for reading normal FASTA formatted text file. """
    def setup(self):
        self.fp = tempfile.NamedTemporaryFile(mode='w+')
        lines = [
            '>seq1 here_is_description1',
            'ATGCATGCATGC',
            '>seq2 here_is_description2',
            'ATGCAAN-ATGCAAA',
            '>marker1 here_is_description_for_marker1',
            'CCCCCCCCCCNNCCC',
            '>marker2 here_is_description_for_marker2',
            '000011110011'
        ]
        self.fp.write('\n'.join(lines))
        self.fp.seek(0)
        self.path = self.fp.name

    def teardown(self):
        self.fp.close()

    def test_fasta_file_to_lists(self):
        """ Tests if fasta_file_to_lists function returns contents of 
        a FASTA file as expected. """
        fasta_d = fasta_file_to_lists(self.path, marker_kw='marker')

        exp_sample = {
            'ids': ['seq1', 'seq2'],
            'descriptions': ['here_is_description1', 'here_is_description2'],
            'sequences': ['ATGCATGCATGC', 'ATGCAAN-ATGCAAA']
        }
        assert fasta_d['sample'] == exp_sample, \
            'Sequence data read from file is not the same as expected: {}'\
                .format(fasta_d['sample'])

        exp_marker = {
            'ids': ['marker1', 'marker2'],
            'descriptions': [
                'here_is_description_for_marker1', 
                'here_is_description_for_marker2'],
            'sequences': ['CCCCCCCCCCNNCCC', '000011110011']
        }
        assert fasta_d['marker'] == exp_marker, \
            'Marker data read from file is not the same as expected: {}'\
                .format(fasta_d['marker'])

    @raises(AssertionError)
    def test_alignment_file_assertion(self):
        """ Tests if alignment_file_to_lists function detects error in sequence 
        lengths. """
        fasta_d = alignment_file_to_lists(self.path, marker_kw='marker')

class TestAlignmentFileReader:
    """ Unit tests for functions for reading alignment file. """
    def setup(self):
        """ Create a dummy alignment file. """
        self.fp = tempfile.NamedTemporaryFile(mode='w+')
        lines = [
            '>seq1 here_is_description1',
            'ATGCATGCATGC',
            '>seq2 here_is_description2',
            'ATGCAAN-ATGC',
            '>marker1 here_is_description_for_marker1',
            'CCCCCCCCCCNN',
            '>marker2 here_is_description_for_marker2',
            '000011110011'
        ]
        self.fp.write('\n'.join(lines))
        self.fp.seek(0)
        self.path = self.fp.name


    def teardown(self):
        self.fp.close()

    def test_alignment_file_to_lists(self):
        """ Tests if alignment_file_to_lists function returns contents of 
        a FASTA file as expected. """
        fasta_d = alignment_file_to_lists(self.path, marker_kw='marker')

        exp_sample = {
            'ids': ['seq1', 'seq2'],
            'descriptions': ['here_is_description1', 'here_is_description2'],
            'sequences': ['ATGCATGCATGC', 'ATGCAAN-ATGC']
        }
        assert fasta_d['sample'] == exp_sample, \
            'Sequence data read from file is not the same as expected: {}'\
                .format(fasta_d['sample'])

        exp_marker = {
            'ids': ['marker1', 'marker2'],
            'descriptions': [
                'here_is_description_for_marker1', 
                'here_is_description_for_marker2'],
            'sequences': ['CCCCCCCCCCNN', '000011110011']
        }
        assert fasta_d['marker'] == exp_marker, \
            'Marker data read from file is not the same as expected: {}'\
                .format(fasta_d['marker'])
