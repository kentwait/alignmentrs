# import os
# from alignmentrs.aln import Alignment
# from alignmentrs.alnset import AlignmentSet


# class TestAlignmentSetGetters:
#     # @classmethod
#     # def setup_class(cls):
#     #     # Create a temporary directory
#     #     cls.temp_dir = 'set'
#     #     os.mkdir(cls.temp_dir)
#     #     # Populate temporary directory
#     #     for i in range(3):
#     #         fname = '{}.fa'.format(i)
#     #         path = os.path.join(cls.temp_dir, fname)
#     #         with open(path, 'w') as f:
#     #             print('>seq1', file=f)
#     #             print('ATG'*30, file=f)
#     #             print('>seq2', file=f)
#     #             print('GAC'*30, file=f)
#     #             print('>seq3', file=f)
#     #             print('CGA'*30, file=f)

#     # @classmethod
#     # def teardown_class(cls):
#     #     # Removes all files in the temp dir
#     #     for fname in os.listdir(cls.temp_dir):
#     #         path = os.path.join(cls.temp_dir, fname)
#     #         os.remove(path)
#     #     # Remove temporary directory
#     #     os.rmdir(cls.temp_dir)

#     def setup(self):
#         self.aln_list = []
#         self.aln_set = AlignmentSet('test_set', self.aln_list, metadata=None)

#     def teardown(self):
#         pass

#     def test_nalns_type(self):
#         pass

#     def test_nalns_value(self):
#         pass

#     def test_alignment_names_type(self):
#         pass

#     def test_alignment_names_value(self):
#         pass

#     def test_alignments_type(self):
#         pass

#     def test_alignments_value(self):
#         pass

#     def test_get_item(self):
#         pass

#     def test_del_item(self):
#         pass

#     def test_repr(self):
#         pass

#     def test_len(self):
#         pass


# class TestAlignmentSetConcatenate:
#     def setup(self):
#         pass

#     def teardown(self):
#         pass

#     def test_new_ids(self):
#         pass

#     def test_new_descriptions(self):
#         pass

#     def test_new_sequences(self):
#         pass

#     def test_new_length(self):
#         pass


# class TestAlignmentSetResample:
#     def setup(self):
#         pass

#     def teardown(self):
#         pass

#     def test_resample_0(self):
#         pass

#     def test_resample_1(self):
#         pass

#     def test_resample_with_replacement(self):
#         pass

#     def test_resample_without_replacement(self):
#         pass
