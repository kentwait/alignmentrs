#! /usr/env/bin/ Python3
from alignmentrs_imports import *
# TODO - check valid ids, check invalid ids, 
                # check mixed type input (list, int, str)
                # check invalid type
                # check invalid optional arg. input
#----------------------------------------------------
class TestAlignment:
    """This is a TestClass for alignmentrs v0.8 library
    
    test type: nosetests
    input file name: test_alignment.txt This is an 
    alignment text file
    
    strategy
    --------
    - created dictionary from test file
        dict.keys = sample ids, dict.values = [sample_description, sample_sequence]
        this enables comparisons between dict items 
        (expected outputs) and alignentrs output
    
    tests applied
    -------------
    tests for expected output for valid input #(EVI)
    tests for output type (EOT)
    tests for expected errors (EE)"""
    
    # basic test for expected outputs
    def setup(self):
        # initiates alignment object for tests
        self.aln_file = Alignment.from_fasta('test_alignment.txt', 'test_align', 'marker')
        self.test_aln_file = fasta_to_dict('test_alignment.txt') # for test_dict
        #for methods that exclude marker sequence
        self.test_aln_file_wo_marker =  fasta_to_dict('test_alignment.txt')
        del self.test_aln_file_wo_marker['marker_0']

    def teardown(self):
        pass
#----------------------------------------------------    
    def test_nrows(self):
        """checks if algn_obj.nrow returns expected type 
        and number of rows in alignment. 
        compares number of keys in test_dict(sample_ids) 
        against algn_obj.nrow output.
        
        pass if:
            algn_obj.nrow returns int and matches 
            len(dict.keys)"""
        
        nrows = self.aln_file.nrows
        test_nrows = len(list(row for row in self.test_aln_file.keys())) #checks number of ids in dict.keys
        
        assert isinstance(nrows,int) #(EOT)
        assert nrows == test_nrows
#----------------------------------------------------    
    def test_nsamples(self):   
        """checks if algn_obj.nsamples returns expected type 
        and number of samples in alignment.
        compares algn_obj.nsamples against len(dict.items).
        
        pass if:
            algn_obj.nrow returns int and matches 
            len(dict.items)"""
        
        nsamples = self.aln_file.nsamples
        test_nsamples = len(list(sample for sample in \
                                 self.test_aln_file_wo_marker.items()))  
        
        assert isinstance(nsamples,int) #(EOT) 
        assert nsamples == test_nsamples
#----------------------------------------------------     
    def test_nmarkers(self):   
        """checks if algn_obj.nmarkers returns expected 
        output type and number of markers in alignment.
        compares number of markers in algn_obj. against 
        number of test_dict keys with 'marker' keyword
        
        pass if:
         output is int
         number of markers in aln.obj == number of markers
         in test_dict"""
        
        nmarkers = self.aln_file.nmarkers
        test_nmarkers = len(list(sample for sample in self.test_aln_file.items() \
                                 if 'marker' in sample[0]))
        
        assert  isinstance(nmarkers, int)
        assert  nmarkers == test_nmarkers
#----------------------------------------------------     
    def test_nsites(self):
        """checks if algn_obj.nsites returns expected 
        output type and number of sites in alignment.
        compares number of sites in algn_obj. against 
        sequence length of a random key in test_dict.
        *note: assumes all samples have same len
        
        pass if:
         output is int
         number of sites in aln.obj == sequence length
         of any test_dict.key()"""
        
        nsites = self.aln_file.nsites
        # extracts random test_dict id to use for seq len
        random_seq_id = random.choice(list(self.test_aln_file.keys()))
        random_seq_len = self.test_aln_file[random_seq_id][1] 
        test_nsites = len(random_seq_len)
        
        assert isinstance(nsites, int)
        assert nsites == test_nsites
#----------------------------------------------------     
    def test_sample_ids(self):
        """checks if all sample ids in align.obj matches
        ids in test_dictr. checks if method returns
        expected type: (list)
        pass if:
            sample_ids in align.obj matches keys in test_dict""" 
        sample_ids = self.aln_file.sample_ids
        test_sample_ids = list(sample_id for sample_id in \
                               self.test_aln_file.keys() if not 'marker' in sample_id)
        
        assert isinstance(sample_ids,list)
        assert sample_ids == test_sample_ids
#----------------------------------------------------     
    def test_sample_descriptions(self):   
        """checks if descriptions in align.obj match descriptions
        in test_dict. checks if expected type is returned
        pass if:
          expected type : list
          descriptions in align.obj matches descriptions 
          in test_dict."""
        
        sample_descriptions =  self.aln_file.sample_descriptions
        # description is first index pos in test_dict.values
        test_descriptions = list(self.test_aln_file[desc][0] for desc in \
                                 self.test_aln_file.keys() if 'marker' not in desc) 
        
        assert  isinstance(sample_descriptions,list)
        assert  sample_descriptions ==  test_descriptions
#----------------------------------------------------     
    def test_sample_sequences(self):   
        """checks if sequences in align.obj match sequences
        in test_dict. checks if expected type is returned
        pass if:
          expected type : list
          sequences in align.obj matches sequences 
          in test_dict."""
        
        sample_sequences =  self.aln_file.sample_sequences
        # sequence is second index pos in test_dict.values
        test_sample_sequences =  list(self.test_aln_file[sample_id][1] for sample_id in \
                                      self.test_aln_file.keys() if not 'marker' in sample_id)
        
        assert isinstance(sample_sequences,list)
        assert sample_sequences == test_sample_sequences
#----------------------------------------------------        
    def test_marker_ids(self):
        """checks if marker_ids in align.obj match marker_ids
        in test_dict. checks if expected type is returned
        pass if:
          expected type : list
          marker_ids in align.obj matches marker_ids 
          in test_dict."""
        
        marker_ids =  self.aln_file.marker_ids
        # marker ids are considered as any key that has 'marker keyword in test_dict'
        test_marker_ids = list(marker_id for marker_id in \
                               self.test_aln_file.keys() \
                               if 'marker' in marker_id)
        
        assert isinstance(marker_ids,list)
        assert marker_ids == test_marker_ids
#----------------------------------------------------       
    def test_marker_descriptions(self):
        """checks if marker_descriptions in align.obj and
        test_dict match. checks if expected type is returned
        pass if:
          expected type : list
          marker_descriptions in align.obj matches marker_descriptions 
          in test_dict."""
        
        marker_descriptions =  self.aln_file.marker_descriptions
        test_marker_descriptions = list(self.test_aln_file[marker_desc][0] \
                                    for marker_desc in self.test_aln_file.keys() \
                                    if 'marker' in marker_desc)
        
        assert isinstance(marker_descriptions,list)
        assert marker_descriptions == test_marker_descriptions
#----------------------------------------------------       
    def test_marker_sequences(self): 
        """checks if marker_sequences in align.obj and
        test_dict match. checks if expected type is returned
        pass if:
          expected type : list
          marker_sequences in align.obj matches marker_sequences 
          in test_dict."""
        
        marker_sequences =  self.aln_file.marker_sequences
        test_marker_sequences = list(self.test_aln_file[marker_id][1] \
                                        for marker_id in self.test_aln_file.keys() \
                                        if 'marker' in marker_id)
        
        assert isinstance(marker_sequences,list)
        assert marker_sequences == test_marker_sequences
#----------------------------------------------------     
    # class methods
    
    # Does not work properly
    
    #def test_subset(cls):
     #   """Returns a subset of the alignment by samples, markers and sites."""
     #   sample_ids = [1,2,3]
     #   marker_ids = ['marker_0']
      #  aln_obj = Alignment.from_fasta('test_alignment.txt', 'test_align', 'marker')
        
      #  subset = (aln_obj,sample_ids, marker_ids)
      #  test_subset = 
#----------------------------------------------------          
    def test_get_samples(self):
        # TODO - check valid ids, check invalid ids, 
                # check mixed type input (list, int, str)
                # check invalid type
                # check invalid optional arg. input
               
        """Returns a list of sequence strings containing only the samples
        specified by the index."""
        test_names = ['Dmel_RG4N', 'Dmel_RG7']
        get_samples = self.aln_file.get_samples(test_names) 
            
        test_list = []
        #extracts expected sequences from test dictionary
        for index in range(len(self.test_aln_file['Dmel_RG4N'][1])):
            nuc1 = self.test_aln_file['Dmel_RG4N'][1][index]
            nuc2 = self.test_aln_file['Dmel_RG7'][1][index]
            test_list = test_list + [[nuc1,nuc2]]
            
        assert [sample for sample in get_samples] == test_list
#----------------------------------------------------     
    # can't test yet becausse returns an object
    #def test_get_markers(self, i, match_prefix=False, match_suffix=False):
            #"""Returns a list of sequence strings containing only the markers
            #specified by the index."""
            
    # error with row and column posisiton
    #def test_get_sites(self):
            #"""Returns a new alignment containing only the sites specified
            #by the given list of column numbers."""
#----------------------------------------------------  
        # Setter/Replacer
    def test_replace_samples(self):
        """Replaces the sequence for a given row in the alignment matrix."""
        new_sample = 'A' * 3432
        replace_sample = self.aln_file.replace_samples(['Dmel_528_2597'], [new_sample])
        sample_sequences =  self.aln_file.sample_sequences
        
        assert new_sample in sample_sequences
#----------------------------------------------------  
    def test_insert_samples_from_lists(self):
        """Inserts new sequences in the alignment matrix at the specifiedrow position inplace."""
    
        new_sample_id = 'Dsim99201'
        new_sample_desc = '|CH2912|'
        new_sample_sequence = 'GCCGATGT' * 429
        self.aln_file.insert_samples_from_lists(1, [new_sample_id], [new_sample_desc], [new_sample_sequence])
        #assert self.aln_file.sample_sequences[1:2] == new_sample_sequence # error here
        assert self.aln_file.sample_ids[1] == new_sample_id 
        assert self.aln_file.sample_descriptions[1] == new_sample_desc 
#----------------------------------------------------  
    def test_append_sample_from_lists(self):
        """Inserts new sequences after the last row of the alignment matrix
        inplace. This increases the total number of samples."""
            
        new_sample_id = 'Dere_lastind'
        new_sample_desc = '|last_ind|'
        new_sample_sequence = 'GCCGAAAA' * 429
        self.aln_file.append_sample_from_lists([new_sample_id], [new_sample_desc], [new_sample_sequence])
        assert self.aln_file.sample_sequences[-1] == new_sample_sequence 
        assert self.aln_file.sample_ids[-1] == new_sample_id 
        assert self.aln_file.sample_descriptions[-1] == new_sample_desc 
#----------------------------------------------------   
    def test_remove_samples(self):
        """Removes sample sequences based on the given index."""
        index_to_remove = 1
        sequence_to_remove = self.aln_file.sample_sequences[index_to_remove]
        self.aln_file.remove_samples(index_to_remove)
        assert  sequence_to_remove != self.aln_file.sample_sequences[1] 


    def test_retain_samples(self):
        """Keeps sample sequences based on the given index."""
        
        self.aln_file.retain_samples([1,2,3,4,5])
        assert self.aln_file.nsamples == 5
        assert len(self.aln_file.sample_ids) == 5
        assert len(self.aln_file.sample_descriptions) == 5
        assert len(self.aln_file.sample_sequences) == 5
        assert self.aln_file.nmarkers == 1
        assert len(self.aln_file.marker_ids) == 1
        assert len(self.aln_file.marker_descriptions) == 1
        assert len(self.aln_file.marker_sequences) == 1
#----------------------------------------------------               
    def test_remove_sites(self):
        """Removes sites based on the given list of column numbers."""
        self.aln_file.remove_sites([0,1,2])
        
        assert self.aln_file.nsites == (3432-3)
        #TODO- should assert to that the positions are not ATG once 
        # get samples sequences is fixes
#----------------------------------------------------   
    def test_retain_sites(self):
        """Keeps sites based on the given list of column numbers."""
        self.aln_file.retain_sites([0,1,2])
        assert self.aln_file.nsites == 3
        #TODO- should assert to that the positions are not ATG once 
        # get samples sequences is fixes
#----------------------------------------------------                
    def test_from_fasta(cls):
        """Create an Alignment object from a FASTA-formatted file."""
        aln_obj = Alignment.from_fasta('test_alignment.txt', 'test_align', 'marker')
        assert aln_obj.nsites == 3432
        assert aln_obj.nsamples == 65
        assert aln_obj.nmarkers == 1
#----------------------------------------------------               
        # Format converters
    # tested, remove hashtag when done
    #def test_to_fasta(self):
        """Saves the alignment as a FASTA-formatted text file."""
        #create_aln_file = self.aln_file.to_fasta('test_aln_file') # should change input to path_to_file
        #TODO- check if len(currentdir) increases by one, 
        #check if test_aln_file in dir
        #use from_fasta to create new alingment object 
#----------------------------------------------------               
    def test_to_sample_matrix(self):
        """Converts sequences into a numpy matrix.""" # is this matrix or array?
        seq_matrix = self.aln_file.to_sample_matrix()
        test_seq_matrix = array([list(key[1]) for key in self.test_aln_file.values() if 'A' in key[1]]) # use key names
        assert  array_equal(seq_matrix,  test_seq_matrix)
#----------------------------------------------------   
    def test_to_marker_matrix(self, size=1):
        """Converts sequences into a numpy matrix."""
        mseq_matrix = self.aln_file.to_marker_matrix()
        test_mseq_matrix = array([list(self.test_aln_file[marker][1]) for marker in \
                                  self.test_aln_file.keys() if 'marker' in marker])
        assert  array_equal(mseq_matrix, test_mseq_matrix)
#----------------------------------------------------   
        # Iterators
    def test_iter_sites(self):
        """Iterates column-wise over the alignment"""
        iter_seq = list(self.aln_file.iter_sites(0,1))
            
        test_cols = []
        seqs = []
        for num in range(0,1):
            for value in self.test_aln_file.values(): 
                seq = value[1]
                #if not 'A' in seq:
                      #continue
                seqs.append(seq[num])

            test_cols.append(sorted(seqs))
            seqs = []
                
        print(test_cols, iter_seq) 
        assert test_cols ==  iter_seq
#----------------------------------------------------   
    def test_iter_sample_sites(self, start=0, stop=None, size=1):
        """Iterates column-wise over the sample alignment. Excludes markers. """
        iter_seq = list(self.aln_file.iter_sample_sites(0,10))

        test_cols = []
        seqs = []
        for num in range(0,10):
            for value in self.test_aln_file_wo_marker.values(): 
                seq = value[1]
                seqs.append(seq[num])

            test_cols.append(seqs)
            seqs = []
               
        assert test_cols ==  iter_seq
#----------------------------------------------------   
    def test_iter_marker_sites(self):
        """Iterates column-wise over the marker alignment. Excludes samples. """
        iter_seq = list(self.aln_file.iter_marker_sites(0,10))

        test_cols = []
        seqs = []
        for num in range(0,10):
            for value in self.test_aln_file.values():
                    
                seq = value[1]
                if set(seq) != {'C', 'N'}:
                    continue
                seqs.append(seq[num])

            test_cols.append(seqs)
            seqs = []
               
        assert test_cols ==  iter_seq
#----------------------------------------------------               
            # How does this work?
 
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