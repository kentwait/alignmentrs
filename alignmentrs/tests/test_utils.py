""" Unit test for utils functions. """

import tempfile

from alignmentrs.utils import fasta_file_to_lists

class TestFastaFileReader:
    """ Unit tests for fasta_file_to_lists """
    def setup(self):
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

    def test_fasta_file_to_lists(self):
        """ Tests if fasta_file_to_lists function returns contents of 
        a FASTA file as expected. """
        fasta_d = fasta_file_to_lists(self.path, marker_kw='marker')

        exp_sample = {
            'ids': ['seq1', 'seq2'],
            'descriptions': ['here_is_description1', 'here_is_description2'],
            'sequences': ['ATGCATGCATGC', 'ATGCAAN-ATGC']
        }
        assert fasta_d['sample'] == exp_sample, \
            'Sequence data read from file is not the same as expected: {}'\
                .foramt(fasta_d['sample'])

        exp_marker = {
            'ids': ['marker1', 'marker2'],
            'descriptions': [
                'here_is_description_for_marker1', 
                'here_is_description_for_marker2'],
            'sequences': ['CCCCCCCCCCNN', '000011110011']
        }
        assert fasta_d['marker'] == exp_marker, \
            'Marker data read from file is not the same as expected: {}'\
                .foramt(fasta_d['marker'])

