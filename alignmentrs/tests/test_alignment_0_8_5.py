import os
from alignmentrs.aln import Alignment

#TODO- aln.subset: test invalid inputs
#      aln.subset: test int and str input  
#      iter__xxx: write tests
#      to_xxx_matrix_: write tests 
#      block_functions: write tests
#      description editing

def type_error(expected, actual):
    return 'Expected type {}, instead got {}'.format(expected, actual)

def value_error(expected, actual):
    return 'Expected value {}, instead got {}'.format(expected, actual)

def index_error(expected, actual):
    return 'Expected type {}, instead got {}'.format(expected, actual)

TypeError_msg = 'Expected TypeError'
ValueError_msg = 'Expected ValueError'
IndexError_msg = 'Expected IndexError'
Error_msg = 'Expected Error'

class TestAlignmentGetters:
    
    # TODO: Refactor tests so that each test is very simple and
    # multiple assertions are minimized, see test_nrows* and test_nsamples*
    # TODO: Make "printer" functions like type_error and value_error
    # above to standardize errors when the assertion
    
    # basic test for expected outputs
    def setup(self):
        # Create an alignment for testing
        self.temp_filename = 'temp.aln'
        with open(self.temp_filename, 'w') as fp:
            print('>marker_0 |91 sp|\n', file=fp)
            print('CCCCCCCCCCCCCCCCCCCCCCCCCC\n', file=fp)
            print('>Dmel_528_2597 |10 sp|\n', file=fp)
            print('ATGAAGAGCAAGGTGGGGGGGGGGGG\n', file=fp)
            print('>Dmel_RG2 |47 sp|\n', file=fp)
            print('ATGAAGAGCAAGGTGGACCCCCCCCC\n', file=fp)
            print('>Dmel_RG4N |15 sp|\n', file=fp)
            print('ATGAAGAGCAAGGTGGAAAAAAAAAA\n', file=fp)
        # initiates alignment object for tests
        self.aln_file = Alignment.from_fasta(self.temp_filename, 'test_align', marker_kw='marker')

    def teardown(self):
        if os.path.exists(self.temp_filename):
            os.remove(self.temp_filename)

    def test_nrows_type(self):
        """Checks if aln.obj.nrows output type matches expected 
        type
        """
        expected = int
        result = type(self.aln_file.nrows)
        assert expected == result, type_error(expected, result)

    def test_nrows_value(self):
        """Checks if aln.obj.nrows output value matches expected 
        value
        """
        expected = 4
        result = self.aln_file.nrows
        assert expected == result, value_error(expected, result)
#-----------------------------------------------------------------------------
    def test_nsamples_type(self):
        """Checks if aln.obj.nrows output type matches expected
       type
        """
        expected = int
        result = type(self.aln_file.nsamples)
        assert expected == result, type_error(expected, result)

    def test_nsamples_value(self):
        """Checks if aln.obj.nrows output value matches expected 
        value
        """
        expected = 3
        result = self.aln_file.nsamples
        assert expected == result, value_error(expected, result)
#-----------------------------------------------------------------------------
    def test_nmarkers_type(self):
        """Checks if aln.obj.nmarkers output type matches
        expected type
        """
        expected = int
        result = self.aln_file.nmarkers
        assert isinstance(result, expected), type_error(expected, result)
        
    def test_nmarkers_value(self):
        """Checks if aln.obj.nmarkers value matches
        expected number of markers in the sample
        """
        expected = 1
        result = self.aln_file.nmarkers
        assert expected == result, value_error(expected, result) 
#-----------------------------------------------------------------------------
    def test_nsites_type(self):
        """checks if aln.obj.nsites output type matches
        expected type
        """
        expected = int
        result = self.aln_file.nsites
        assert isinstance(result, expected), type_error(expected, result)
        
    def test_nsites_value(self):
        """checks if aln.obj.nsites value matches expected number
        of sites in alignment(i.e sequence length) 
        """
        expected = 26
        result = self.aln_file.nsites
        assert expected == result, value_error(expected, result)
#-----------------------------------------------------------------------------        
    def test_sample_ids_type(self):
        """Checks if aln.obj.sample_ids type matches
        expected type
        """
        expected = list
        result = self.aln_file.sample_ids
        assert isinstance(result, expected), type_error(expected, result) 
    
    def test_sample_ids_value(self):
        """Checks if all aln.obj.sample_ids match expected 
        sample ids
        """
        expected = ['Dmel_528_2597', 'Dmel_RG2', 'Dmel_RG4N']
        result = self.aln_file.sample_ids
        assert expected == result, value_error(expected, result)  
#-----------------------------------------------------------------------------    
    def test_sample_descriptions_type(self):
        """Checks if aln.obj.sample_descriptions output type 
        matches expected type
        """
        expected = list
        result = self.aln_file.sample_descriptions
        assert  isinstance(result, expected), type_error(expected, result) 
        
    def test_sample_descriptions_value(self):
        """Checks if aln.obj.sample_descriptions value matches 
        expected value
        """
        expected = ['|10 sp|', '|47 sp|', '|15 sp|']
        result = self.aln_file.sample_descriptions
        assert  expected == result, value_error(expected, result) 
#-----------------------------------------------------------------------------
    def test_sample_sequences_type(self):
        """Checks if aln.obj.sample_sequences output type 
        matches expected type
        """
        expected = list
        result = self.aln_file.sample_sequences
        assert isinstance(result, expected), type_error(expected, result) 
        
    def test_sample_sequences_value(self):
        """Checks if aln.obj.sample_sequences value output matches 
        expected value
        """
        expected = ['ATGAAGAGCAAGGTGGGGGGGGGGGG',
                                     'ATGAAGAGCAAGGTGGACCCCCCCCC', 
                                     'ATGAAGAGCAAGGTGGAAAAAAAAAA']
        result = self.aln_file.sample_sequences
        assert result == expected, value_error(expected, result)   
#-----------------------------------------------------------------------------           
    def test_marker_ids_type(self):
        """Checks if aln.obj.marker_ids output type matches expected 
        type
        """
        expected = list
        result = self.aln_file.marker_ids
        assert isinstance(result, expected), type_error(expected, result) 
        
    def test_marker_ids_value(self):
        """Checks if aln.obj.marker_ids value matches 
        expected value 
        """
        expected = ['marker_0']
        result = self.aln_file.marker_ids
        assert result == expected, value_error(expected, result)
#-----------------------------------------------------------------------------        
    def test_marker_descriptions_type(self):
        """checks if aln.obj.marker_descriptions output type
        matches expected output
        """ 
        expected = list
        result = self.aln_file.marker_descriptions
        assert isinstance(result, expected), type_error(expected, result)
        
    def test_marker_descriptions_value(self):
        """checks if aln.obj.marker_descriptions output
        value expected output
        """
        expected = ['|91 sp|']
        result = self.aln_file.marker_descriptions
        assert result == expected, value_error(expected, result)
#-----------------------------------------------------------------------------   
    def test_marker_sequences_type(self): 
        """Checks if align.obj.marker_sequences output type
        matches expected output
        """
        expected = list
        result = self.aln_file.marker_sequences
        assert isinstance(result, expected), type_error(expected, result)
        
    def test_marker_sequences_value(self): 
        """Checks if align.obj.marker_sequences value
        matches expected output
        """
        expected = ['CCCCCCCCCCCCCCCCCCCCCCCCCC']
        result = self.aln_file.marker_sequences
        assert result == expected, value_error(expected, result)
#----------------------------------------------------     
    # class methods
    def test_subset_input_sample_id(self):
        """Returns a subset of the alignment by samples, markers and sites."""
        sample_ids = [0,1]
        expected_nsamples = 2
        expected_names = ['Dmel_528_2597', 'Dmel_RG2']
        expected_descriptions = ['|10 sp| ', '|47 sp| ']
        expected_sequences = ['ATGAAGAGCAAGGTGGACCCCCCCCC', 'ATGAAGAGCAAGGTGGAAAAAAAAAA']
        result_subset = self.aln_file.subset(self.aln_file, sample_ids) # fails for sume reason
        result_nsamples = result_subset.nsamples
        result_ids = result_subset.sample_ids
        result_desc = result_subset.sample_descriptions
        result_seq = result_subset.sample_sequences
        assert expected_nsamples == result_nsamples, value_error(expected_nsamples,result_nsamples)
        assert expected_names == result_ids, value_error(expected_names ,result_ids)
        assert expected_descriptions == result_desc, value_error(expected_descriptions ,result_desc)
        assert expected_sequences == result_seq, value_error(expected_sequences ,result_seq)
        
    def test_subset_input_marker_id(self):
        """Returns a subset of the alignment by samples, markers and sites."""
        marker_index = [0]
        expected_nsamples = 1
        expected_names = ['marker_0']
        expected_descriptions = ['|91 sp|']
        expected_sequences = ['CCCCCCCCCCCCCCCCCCCCCCCCCC']
        result_subset = self.aln_file.subset(self.aln_file, marker_ids = marker_index)
        result_nmarkers = result_subset.nmarkers
        result_ids = result_subset.marker_ids
        result_desc = result_subset.marker_descriptions
        result_seq = result_subset.marker_sequences
        assert expected_nsamples == result_nmarkers, value_error(expected_nsamples,result_nmarkers)
        assert expected_names == result_ids, value_error(expected_names ,result_ids) # fails, gives wrong expected marker
        assert expected_descriptions == result_desc, value_error(expected_descriptions ,result_desc)
        assert expected_sequences == result_seq, value_error(expected_sequences ,result_seq)
     
   # HSD: TODO- test invalid invalid inputs 
   # HSD: TODO- test int and str input
    
   # def test_subset(self):
   #     """Returns a subset of the alignment by samples, markers and sites."""
   #     expected_sample_ids = [0,1]
   #     result_subset = (self.aln_file,sample_ids)
        
   #     result_sample_ids = result_subset.get_markers
   #     result_sample_desc = 
   #    result_sample_seq = 
   #     assert
#----------------------------------------------------          
    def test_get_samples_valid_str(self):
        """checks if aln.object.get_samples returns expected
        output type given an input str or int"""
        
        result = self.aln_file.get_samples('Dmel_528_2597')
        result_0 = self.aln_file.get_samples(0)
        expected = object
        assert isinstance(result,expected), type_error(expected, result)
        assert isinstance(result_0,expected), type_error(expected, result_0)
        
    def test_get_samples_valid_list(self):
        """checks if aln.object.get_samples returns expected
        output type given an input list_str or list_int"""
        
        result = self.aln_file.get_samples(['Dmel_528_2597'])
        result_0 = self.aln_file.get_samples([0])
        expected = object
        assert isinstance(result,expected), type_error(expected, result)
        assert isinstance(result_0,expected), type_error(expected, result_0)
#-----------------------------------------------------------------------------        
    def test_get_samples_str_in_invalid_type(self):
        """checks if aln.object.get_samples raises error, given
        an invalid type i.e set, tuple or none"""
        try:
            valid_str_invalid_type_set_str = self.aln_file.get_samples({})
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            valid_str_invalid_type_tuple_str = self.aln_file.get_samples(())
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            None_type = self.aln_file.get_samples(None)
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)    
            
    def test_get_samples_invalid_list_items(self):
        """checks if aln.object.get_samples raises error, given
        an invalid type in a list i.e tuple, set none in list"""
        try:
            set_in_list = self.aln_file.get_samples([{}])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            tuple_in_list = self.aln_file.get_samples([()])
        except TypeError:
            pass      
        else:
            raise Exception(TypeError_msg)
        try:
            none_in_list = self.aln_file.get_samples([None])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
            
    def test_get_samples_invalid_str_sample_id(self):
        """checks if aln.object.get_samples raises error, given
        an invalid string sample id"""
        try:
            invalid_str = self.aln_file.get_samples('Dmel_528_259700')
        except ValueError:
            pass
        else:
            raise Exception(ValueError_msg)
    
    def test_get_samples_invalid_int_sample_id(self):
        """checks if aln.object.get_samples raises error, given
        an invalid index of a sample id"""
        try:
            invalid_str = self.aln_file.get_samples(-1)
        except ValueError:
            pass
        else:
            raise Exception(ValueError_msg)        
            
    def test_get_samples_invalid_str_sample_id(self):
        """checks if aln.object.get_samples raises error, given
        an invalid string sample id in list"""
        try:
            invalid_str = self.aln_file.get_samples(['Dmel_528_259700'])
        except ValueError:
            pass
        else:
            raise Exception(ValueError_msg)
    
    def test_get_samples_invalid_int_sample_id(self):
        """checks if aln.object.get_samples raises error, given
        an invalid index of a sample id in list"""
        try:
            invalid_index = self.aln_file.get_samples([-1])
        except IndexError:
            pass
        else:
            raise Exception(IndexError_msg)   
            
        # TODO: Do not make compound tests!
        # Each assertion should be a separate test so that it will fail for 
        # that specific test.
        # For example, test.nsites == expected_sites is unnecessary because it
        # duplicates an existing test - test_nsites
#-----------------------------------------------------------------------------  
    def test_get_markers_valid_input(self):
            """checks if aln.object.get_markers returns expected
            marker sequences of an alignment file
            """
            result = self.aln_file.get_markers('marker_0')
            result_0 = self.aln_file.get_markers(0)
            expected = object
            assert isinstance(result,expected), type_error(expected, result)
            assert isinstance(result_0,expected), type_error(expected, result_0)
            
    def test_get_markers_valid_list_input(self):
            """Tests if aln.object.get_markers returns expected
            marker sequences of an alignment file, given a list
            input
            """
            result = self.aln_file.get_markers(['marker_0'])
            result_0 = self.aln_file.get_markers([0])
            expected = object 
            assert isinstance(result,expected), type_error(expected, result)
            assert isinstance(result_0,expected), type_error(expected, result_0)
            # TODO: testing single char and list should be different tests as they are two scenarios of calling the method
                      
    def test_get_markers_valid_invalid_type(self):
        """checks if aln.object.get_markers raises error, given
        an invalid type i.e tuple, set none"""
        try:
            valid_str_invalid_type_0 = self.aln_file.get_markers({})
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            valid_str_invalid_type_00 = self.aln_file.get_markers(())
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
            
        try:
            None_type = self.aln_file.get_markers(None)
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)    
    def test_get_markers_valid_invalid_list_items(self):
        """checks if aln.object.get_markers raises error, given
        an invalid type in a list i.e tuple, set or none"""
        try:
            valid_str_invalid_type_0 = self.aln_file.get_markers([{}])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            valid_str_invalid_type_00 = self.aln_file.get_markers([()])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            None_type = self.aln_file.get_markers([None])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
#----------------------------------------------------     
    #def test_get_sites(self): # row error
           # """Returns a new alignment containing only the sites specified
           # #by the given list of column numbers."""
            
#----------------------------------------------------  
        # Setter/Replacer
    #TODO- could test invalid input types 
    def test_insert_samples_from_lists_valid_inputs(self):
        """checks if aln.object.insert_samples_from_lists inserts a 
        new sample into the aln.object
        """
        expected_id = 'Dsim99201'
        expected_desc = '|CH2912|'
        expected_sequence = 'G' * 26
        self.aln_file.insert_samples_from_lists(1, [expected_id],[expected_desc],[expected_sequence])
        result_id = self.aln_file.sample_ids[1]
        result_description = self.aln_file.sample_descriptions[1]
        result_sequence = self.aln_file.sample_sequences[1]
        assert expected_id == result_id, value_error(expected_id,result_id)  
        assert  expected_desc == result_description , value_error(expected_desc,result_description)
        assert result_sequence == expected_sequence, value_error(expected_sequence,result_sequence)
        
    def test_insert_samples_from_lists_valid_invalid_index(self):
        """checks if aln.object.insert_samples_from_lists raises
        index error given an invalid index
        """
        invalid_index = -1
        valid_id = 'Dsim99201'
        valid_desc = '|CH2912|'
        valid_sequence = 'G' * 26
        try:
            self.aln_file.insert_samples_from_lists(invalid_index, [valid_id],[valid_desc], [valid_sequence])
        except IndexError:
            pass
        else:
            raise Exception(IndexError_msg)
        
    def test_insert_samples_from_lists_valid_invalid_id(self):
        """checks if aln.object.insert_samples_from_lists raises
        Type error given an invalid type for sample id
        """
        valid_index = 1
        invalid_id = 1
        valid_desc = '|CH2912|'
        valid_sequence = 'G' * 26
        try:
            self.aln_file.insert_samples_from_lists(valid_index,[invalid_id],[valid_desc],[valid_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        
    def test_insert_samples_from_lists_valid_invalid_description(self):
        """checks if aln.object.insert_samples_from_lists raises
        Type error given an invalid type for sample description
        """
        valid_index = 1
        valid_id = 'Dsim99201'
        invalid_desc = 1
        valid_sequence = 'G' * 26
        try:
            self.aln_file.insert_samples_from_lists(valid_index,[valid_id],[invalid_desc],[valid_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        
    def test_insert_samples_from_lists_valid_invalid_sequence(self):
        """checks if aln.object.insert_samples_from_lists raises
        Type error given an invalid type for sample sequence
        """
        valid_index = 1
        valid_id = 'Dsim99201'
        valid_desc = 1
        invalid_sequence =  26
        try:
            self.aln_file.insert_samples_from_lists(valid_index, [valid_id],[valid_desc],[invalid_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
#-----------------------------------------------------------------------------
    def test_append_sample_from_lists_valid_inputs(self):
        """Tests if append_sample_from_lists appends a
        to the last index of aln.object
        """
        expected_id = 'Dere_lastind'
        expected_desc = '|last_ind|'
        expected_sequence = 'G' * 26
        self.aln_file.append_sample_from_lists([expected_id],[expected_desc], [expected_sequence])
        result_id =  self.aln_file.sample_ids[-1] 
        result_desc = self.aln_file.sample_descriptions[-1]
        result_sequence = self.aln_file.sample_sequences[-1]
        assert result_id == expected_id,  type_error(expected_id, result_id)  
        assert result_desc == expected_desc,  type_error(expected_desc, result_desc)  
        assert result_sequence == expected_sequence,  type_error(expected_sequence, result_sequence)
        
    def test_append_sample_from_lists_invalid_id_type(self):
        """checks if aln.object.append_sample_from_lists raises
        Type error given an invalid type for sample id
        """
        sample_id = 1
        sample_desc = '|last_ind|'
        sample_sequence = 'G' * 26
        
        try:
            result = self.aln_file.append_sample_from_lists([sample_id],[sample_desc], [sample_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        
    def test_append_sample_from_lists_invalid_description_type(self):
        """checks if aln.object.append_sample_from_lists raises
        Type error given an invalid type for sample description
        """
        expected_id = 'Dere_lastind'
        expected_desc = 1
        expected_sequence = 'G' * 26
        try:
            result = self.aln_file.append_sample_from_lists([expected_id],[expected_desc], [expected_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)

    def test_append_sample_from_lists_invalid_sequence_type(self):
        """checks if aln.object.append_sample_from_lists raises
        Type error given an invalid type for sample sequence
        """
        expected_id = 'Dere_lastind'
        expected_desc = '|last_ind|'
        expected_sequence =  26
        try:
            result = self.aln_file.append_sample_from_lists([expected_id],[expected_desc], [expected_sequence])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
#-----------------------------------------------------------------------------            
    def test_remove_samples_valid_sample_name(self):
        """checks if aln.obj.remove_samples removes sample
        information  of a sample in the alignment object
        given a str sample id
        """
        sample_name_to_remove = 'Dmel_RG2'
        sample_index = 1
        self.aln_file.remove_samples(sample_name_to_remove)
        result_sample_id =  self.aln_file.sample_ids[sample_index]
        result_sequence = self.aln_file.sample_sequences[sample_index]
        result_desc  =  self.aln_file.sample_descriptions[sample_index]
        expected_sample_id = ['Dmel_RG4N']
        expected_sequence = ['|15 sp|']
        expected_desc = 'ATGAAGAGCAAGGTGGAAAAAAAAAA'
        assert  expected_sample_id != result_sample_id, value_error(expected_sample_id,result_sequence)
        assert  expected_sequence != result_desc, value_error(expected_sequence,result_sample_id) 
        assert  expected_desc != result_desc,type_error(expected_desc,result_desc)
        
    def test_remove_samples_valid_index(self):
        """checks if aln.obj.remove_samples removes sample
        information  of a sample in the alignment object
        given the sample index position
        """
        index_to_remove = 1
        self.aln_file.remove_samples(index_to_remove)
        result_sample_id =  self.aln_file.sample_ids[index_to_remove]
        result_sequence = self.aln_file.sample_sequences[index_to_remove]
        result_desc  =  self.aln_file.sample_descriptions[index_to_remove]
        expected_sample_id = ['Dmel_RG4N']
        expected_sequence = ['|15 sp|']
        expected_desc = 'ATGAAGAGCAAGGTGGAAAAAAAAAA'
        assert  expected_sample_id != result_sample_id, value_error(expected_sample_id,result_sequence)
        assert  expected_sequence != result_desc, value_error(expected_sequence,result_sample_id ) 
        assert  expected_desc != result_desc,type_error(expected_desc,result_desc ) 
    
    def test_remove_samples_invalid_str(self):
        """checks if aln.obj.remove_samples gives a 
        TypeError given an invalid sample id
        """
        sample_name_to_remove = 'Dmel_RG4N0'
        try:
            self.aln_file.remove_samples(sample_name_to_remove)
        except ValueError:
            pass
        else:
            raise Exception(TypeError_msg) 
            
    def test_remove_samples_invalid_index(self):
        """checks if aln.obj.remove_samples gives an 
        IndexError given an invalid input index
        """
        index_to_remove = -1
        try:
            self.aln_file.remove_samples(index_to_remove)
        except IndexError:
            pass
        else:
            raise Exception(index_error)
   
    def test_remove_samples_invalid_type(self):
        """checks if aln.obj.remove_samples gives a 
        ValueError given an invalid input type i.e 
        dict, float, set and tuple
        """    
        float_type = 1.1
        dict_type = {1:'w'}
        set_type = {}
        tuple_type = ()
        try:
            self.aln_file.remove_samples(float_type) 
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        
        try:
            self.aln_file.remove_samples(dict_type) 
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.remove_samples(set_type) 
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.remove_samples(tuple_type) 
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
#-----------------------------------------------------------------------------            
    def test_retain_samples_valid_input(self):
        """checks if aln.obj.retain_samples changes 
        properties of the alignment object
        """
        self.aln_file.retain_samples([1])
        expected_len = 1
        result_sample_num =  self.aln_file.nsamples 
        result_sample_id_num =  len(self.aln_file.sample_ids)
        result_sample_desc_num =  self.aln_file.sample_descriptions
        result_sample_seq_num = len(self.aln_file.sample_sequences)
        result_marker_num =  self.aln_file.nmarkers
        result_marker_id_num = len(self.aln_file.marker_ids)
        result_marker_desc_num =  len(self.aln_file.marker_descriptions)
        result_marker_seq_num =  len(self.aln_file.marker_sequences)
        
        assert result_sample_num == expected_len, value_error(expected_len, result_sample_num )
        assert result_sample_id_num == expected_len, value_error(expected_len, result_sample_id_num )
        assert result_sample_desc_num == expected_len, valu_error(expected_len, result_sample_desc_num )
        assert result_sample_seq_num  == expected_len, value_error(expected_len, result_sample_seq_num)
        assert result_marker_num == expected_len, value_error(expected_len, result_marker_num ) 
        assert result_marker_id_num == expected_len, value_error(expected_len, result_marker_id_num )
        assert result_marker_desc_num == expected_len, value_error(expected_len, result_marker_desc_num)
        assert result_marker_seq_num == expected_len, value_error(expected_len, result_marker_seq_num)
        
    def test_retain_samples_valid_input(self):
        """checks if aln.obj.retain_samples retains
        expected sample in the alignment
        """    
        self.aln_file.retain_samples([1])
        expected_len = 1
        result_sample_id = self.aln_file.sample_ids[0] 
        result_sample_description = self.aln_file.sample_descriptions[0]
        result_sequence = self.aln_file.sample_sequences[0]
        expected_sample_id =  'Dmel_RG2'
        expected_sample_description = '|47 sp|' 
        expected_sequence = 'ATGAAGAGCAAGGTGGACCCCCCCCC'
        
        assert expected_sample_id ==  result_sample_id,  value_error(expected_sample_id, result_sample_id)
        assert expected_sample_description == result_sample_description,  value_error(expected_sample_description, result_sample_description)
        assert expected_sequence == result_sequence,  value_error(expected_sequence, result_sequence)
 
        def test_retain_samples_invalid_list_input(self):
            """checks if aln.obj.retain_samples raises
            an error given the wrong input type in list
            """
            try:
                self.aln_file.retain_samples(['1.1'])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)    
            try:
                self.aln_file.retain_samples([1.1])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)     
            try:
                self.aln_file.retain_samples([-1])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
            try:
                self.aln_file.retain_samples([{}])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
            try:
                self.aln_file.retain_samples([])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)       
            try:
                self.aln_file.retain_samples([None])
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
                
        def test_retain_samples_invalid_input(self):
            """checks if aln.obj.retain_samples raises
            an error given the wrong input type
            """
            try:
                self.aln_file.retain_samples('1.1')
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)    
            try:
                self.aln_file.retain_samples(1.1)
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)     
            try:
                self.aln_file.retain_samples(-1)
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
            try:
                self.aln_file.retain_samples({})
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
            try:
                self.aln_file.retain_samples()
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)       
            try:
                self.aln_file.retain_samples(None)
            except TypeError:
                pass
            else:
                raise Exception(TypeError_msg)
        # TODO: Other than count, test whether the correct sample was retained
#-----------------------------------------------------------------------------
    def test_remove_sites_list_input(self):
        """tests if aln.obj.remove_sites makes expected change
        to number of sites in alignment given int in list"""
        self.aln_file.remove_sites([2])
        expected = 25
        result = self.aln_file.nsites
        assert expected == result, value_error(expected, actual)
        
    def test_remove_sites_int_input(self):
        """tests if aln.obj.remove_sites makes expected change
        to number of sites in alignment given an input int"""
        self.aln_file.remove_sites(2)
        expected = 25
        result = self.aln_file.nsites
        assert expected == result, value_error(expected, actual)
        
    def test_remove_sites_empty_input(self):
        """tests if aln.obj.remove_sites makes expected change
        to number of sites in alignment given an input int"""
        self.aln_file.remove_sites(())
        expected = 26
        result = self.aln_file.nsites
        assert expected == result, value_error(expected, actual)
        
    def test_remove_sites_empty_list_input(self):
        """tests if aln.obj.remove_sites makes expected change
        to number of sites in alignment given an input int"""
        self.aln_file.remove_sites([])
        expected = 26
        result = self.aln_file.nsites
        assert expected == result, value_error(expected, actual)
        
    def test_remove_sites_invalid_input(self):
        """tests if aln.obj.remove_sites raises an 
        error given invalid input"""
        try:
            self.aln_file.remove_sites(-2)
        except IndexError:
            pass
        else:
            raise Exception(IndexError_msg)
        try:
            self.aln_file.remove_sites({})
        except TypeError:
            pass   
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.remove_sites(None)
        except TypeError:
            pass  
        else:
            raise Exception(TypeError_msg)
        
    def test_remove_sites_invalid_list_input(self):
        """tests if aln.obj.remove_sites raises an 
        error given an invalid list input"""
        try:
            self.aln_file.remove_sites([-2])
        except IndexError:
            pass
        else:
            raise Exception(IndexError_msg)
            
        try:
            self.aln_file.remove_sites([()])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
            
        try:
            self.aln_file.remove_sites([{}])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.remove_sites([None])
        except TypeError:
            pass
        else:
            raise Exception(TypeError_msg)
#-----------------------------------------------------------------------------            
    def test_retain_sites_valid_input(self):
        """tests if aln.obj.retain_sites removes all 
        but specified samples in the alignment object"""
        self.aln_file.retain_sites(0,1,2)
        expected = 26
        result = self.aln_file.nsites
        
        assert expected == result, value_error(expected, result) 
        
    def test_retain_sites_valid_list_input(self):
        """tests if aln.obj.retain_sites removes all 
        but specified samples in the alignment object"""
        self.aln_file.retain_sites([0,1,2])
        expected = 3
        result = self.aln_file.nsites
        
        assert expected == result, value_error(expected, result) 
        
    def test_retain_sites_invalid_input(self):
        """tests if aln.obj.retain_sites all sequences
        apart from specified sequences in alignment sequences"""
        try:
            self.aln_file.retain_sites(-1)
        except IndexError: # Expected index error 
            pass
        else:
            raise Exception(IndexError_msg)
        try:
            self.aln_file.retain_sites((1,))
        except TypeError: # Expected index error 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites({})
        except TypeError: # Expected index error 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites(-1)
        except TypeError: # Expected index error 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites('')
        except TypeError: # Expected index error 
            pass
        else:
            raise Exception(TypeError_msg)
        
    def test_retain_sites_invalid_input(self):
        """tests if aln.obj.retain_sites all sequences
        apart from specified sequences in alignment sequences
        given an input list"""
        try:
            self.aln_file.retain_sites([-1])
        except IndexError: # Expected index error 
            pass
        else:
            raise Exception(IndexError_msg)
        try:
            self.aln_file.retain_sites([(1,)])
        except TypeError: 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites([{}])
        except TypeError: 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites([-1])
        except TypeError: 
            pass
        else:
            raise Exception(TypeError_msg)
        try:
            self.aln_file.retain_sites([''])
        except TypeError: 
            pass
        else:
            raise Exception(TypeError_msg)
#-----------------------------------------------------------------------------            
    def test_from_fasta(cls):
        """Create an Alignment object from a FASTA-formatted file."""
        aln_obj = Alignment.from_fasta('temp.aln', 'test_align_from_faster', 'marker')
        result_nsites = aln_obj.nsites
        result_samples = aln_obj.nsamples
        results_nmarkers = aln_obj.nmarkers
        expected_nsites = 26
        expected_nsamples = 3
        expected_nmarkers = 1
        assert result_nsites == expected_nsites,value_error(expected_nsites, result) 
        assert result_samples == expected_nsamples,value_error(expected_nsamples, result)  
        assert results_nmarkers == expected_nmarkers,value_error(expected_nmarkers, result)
#-----------------------------------------------------------------------------
        # Format converters
    # tested, remove hashtag when done
    #def test_to_fasta(self):
#        """Saves the alignment as a FASTA-formatted text file."""
        #create_aln_file = self.aln_file.to_fasta('test_aln_file') # should change input to path_to_file
        #TODO- check if len(currentdir) increases by one, 
        #check if test_aln_file in dir
        #use from_fasta to create new alingment object 
#-----------------------------------------------------------------------------     
#----------------------------------------------------               
#    def test_to_sample_matrix(self):
#        """Converts sequences into a numpy matrix.""" # is this matrix or array?
#        seq_matrix = self.aln_file.to_sample_matrix()
#        test_seq_matrix = array([list(key[1]) for key in self.test_aln_file.values() if 'A' in key[1]]) # use key names
#        assert  array_equal(seq_matrix,  test_seq_matrix)
#----------------------------------------------------   
#    def test_to_marker_matrix(self, size=1):
#        """Converts sequences into a numpy matrix."""
#        mseq_matrix = self.aln_file.to_marker_matrix()
#        test_mseq_matrix = array([list(self.test_aln_file[marker][1]) for marker in \
#                                  self.test_aln_file.keys() if 'marker' in marker])
#        assert  array_equal(mseq_matrix, test_mseq_matrix)
#----------------------------------------------------   
#        # Iterators
#    def test_iter_sites(self):
#        """Iterates column-wise over the alignment"""
#        iter_seq = list(self.aln_file.iter_sites(0,1))
#            
#        test_cols = []
#        seqs = []
#        for num in range(0,1):
#            for value in self.test_aln_file.values(): 
#                seq = value[1]
#                #if not 'A' in seq:
#                      #continue
#                seqs.append(seq[num])
#
#            test_cols.append(sorted(seqs))
#            seqs = []
#                
#        print(test_cols, iter_seq) 
#        assert test_cols ==  iter_seq
#----------------------------------------------------   
##    def test_iter_sample_sites(self, start=0, stop=None, size=1):
##        """Iterates column-wise over the sample alignment. Excludes markers. """
##        iter_seq = list(self.aln_file.iter_sample_sites(0,10))
##
##        test_cols = []
##        seqs = []
##        for num in range(0,10):
##            for value in self.test_aln_file_wo_marker.values(): 
##                seq = value[1]
##                seqs.append(seq[num])
##
##            test_cols.append(seqs)
##            seqs = []
##               
##        assert test_cols ==  iter_seq
#----------------------------------------------------   
###    def test_iter_marker_sites(self):
###        """Iterates column-wise over the marker alignment. Excludes samples. """
###        iter_seq = list(self.aln_file.iter_marker_sites(0,10))
###
###        test_cols = []
###        seqs = []
###        for num in range(0,10):
###            for value in self.test_aln_file.values():
###                    
###                seq = value[1]
###                if set(seq) != {'C', 'N'}:
###                    continue
###                seqs.append(seq[num])
###
###            test_cols.append(seqs)
###            seqs = []
###               
###        assert test_cols ==  iter_seq
#----------------------------------------------------    
# How does this work? Not yet available for v.85
 
        # Block-related methods

   #def test_set_blocklists(self, ref_seq, description_encoder=None):
#            """Creates new block information for the sequences given a reference. """
#----------------------------------------------------   
   # def test_parse_description_as_blocks(self, description_decoder=None):
#            """Parses sample description into block data."""
#----------------------------------------------------   
   # def test_write_blocks_to_description(self, description_encoder):
#            """Writes each sample's block data as a string, replacing its
 #           description."""
#----------------------------------------------------   