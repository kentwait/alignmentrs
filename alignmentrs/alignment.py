from libalignmentrs.alignment import BaseAlignment
import blockrs


class Alignment:
    """Alignment is a complete representation of a multiple sequence alignment
    of biological sequences an their annotations such as alignment markers and
    alignment block data.
    """
    def __init__(self, name, sample_alignment, marker_alignment):
        """Creates a new Alignment object from sample and marker alignments.

        Parameters
        ----------
        sample_alignment
        marker_alignment
        """
        self.name = name
        self.samples = sample_alignment
        self.markers = marker_alignment
        self.blocklists = []
        if not (self.markers is None or self.markers.nsamples == 0):
            assert self.samples.nsites == self.markers.nsites, \
                "Sample and marker nsites are not equal."

    @property
    def nsites(self):
        """Returns the number of sites in the alignment.
        """
        return self.samples.nsites

    @property
    def nsamples(self):
        """Returns the number of samples in the alignment.
        """
        return self.samples.nsamples

    @property
    def nmarkers(self):
        """Returns the number of markers in the alignment
        """
        if self.markers is None:
            return None
        return self.markers.nsamples

    # @property
    # def sample_matrix(self):
    #     """Returns the sample alignent matrix
    #     """
    #     return self._sample_aln.matrix

    # @property
    # def marker_matrix(self):
    #     """Returns the marker alignment matrix
    #     """
    #     return self._marker_aln.matrix

    # @classmethod
    # def subset(cls, aln, sample_ids=None, marker_ids=None, sites=None,
    #            sample_id_step=1, marker_id_step=1, site_step=1):
    #     """Returns a subset of the alignment by samples, markers and sites.

    #     Parameters
    #     ----------
    #     aln : Alignment
    #     sample_ids : list
    #     marker_ids : list
    #     sites : list
    #     sample_id_step : int
    #     marker_id_step : int
    #     site_step : int

    #     Returns
    #     -------
    #     Alignment

    #     """
    #     if sample_ids is None:
    #         sample_ids = range(0, aln.nsamples, sample_id_step)
    #     else:
    #         if sample_id_step != 1:
    #             raise ValueError('sample_id_step value is considered only ' \
    #                              'if sample_ids is None')
    #     if marker_ids is None:
    #         marker_ids = range(0, aln.nmarkers, marker_id_step)
    #     else:
    #         if marker_id_step != 1:
    #             raise ValueError('marker_id_step value is considered only ' \
    #                              'if marker_ids is None')
    #     if sites is None:
    #         sites = range(0, aln.nsites, site_step)
    #     else:
    #         if site_step != 1:
    #             raise ValueError('site_step value is considered only ' \
    #                              'if sites is None')
    #     new_aln = cls.__new__(cls)
    #     new_aln.name = aln.name
    #     new_aln._sample_aln = aln._sample_aln.__class__.subset(
    #         aln._sample_aln,
    #         rows=sample_ids, cols=sites,
    #         row_step=sample_id_step, col_step=site_step
    #     )
    #     new_aln._marker_aln = aln._marker_aln.__class__.subset(
    #         aln._marker_aln,
    #         rows=marker_ids, cols=sites,
    #         row_step=marker_id_step, col_step=site_step
    #     )
    #     return new_aln

    def replace_samples(self, i, sequences):
        """Replaces the sequence for a given row in the alignment matrix.

        Parameters
        ----------
        i : int, str, list
        sequences : str

        """
        if isinstance(i, int) and isinstance(sequences, str):
            self.samples.set_sequences([i], [sequences])
        elif isinstance(i, str) and isinstance(sequences, str):
            ids = self.samples.sample_names_to_ids([i])
            self.samples.set_sequences([ids], [sequences])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            self.samples.set_sequences(i, sequences)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            ids = self.samples.sample_names_to_ids(i)
            self.samples.set_sequences(ids, sequences)
        else:
            raise ValueError('i must be an int, str, list of int, or list of str.')

    # def insert_samples(self, sequence_str, i):
    #     """Inserts a new sequence in the alignment matrix at the specified
    #     row position. This increases the total number of rows.

    #     Parameters
    #     ----------
    #     sequence_str : str or list of str
    #     i : int or list of int

    #     """
    #     self._sample_aln.insert_samples(sequence_str, i)

    # def append_sample(self, sequence_str):
    #     """Inserts a new sequence after the last row of the alignment matrix.
    #     This increases the total number of rows by 1.

    #     Parameters
    #     ----------
    #     sequence_str : str

    #     """
    #     self._sample_aln.append_sample(sequence_str)

    def remove_samples(self, i, match_prefix=False, match_suffix=False):
        """Removes sample sequences based on the given index.
        If index is a number, only one sequence is removed.
        If the index is a list of numbers, the sequence found at each row
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        if isinstance(i, int):
            self.samples.remove_samples([i])
        elif isinstance(i, str):
            if match_prefix:
                self.samples.remove_samples_by_prefix([i])
            elif match_suffix:
                self.samples.remove_samples_by_suffix([i])
            else:
                self.samples.remove_samples_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            self.samples.remove_samples(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                self.samples.remove_samples_by_prefix(i)
            elif match_suffix:
                self.samples.remove_samples_by_suffix(i)
            else:
                self.samples.remove_samples_by_name(i)        
        else:
            raise ValueError('i must be an int, str, list of int, or list of str.')

    def retain_samples(self, i, match_prefix=False, match_suffix=False):
        """Keeps sample sequences based on the given index.
        """
        if isinstance(i, int):
            self.samples.retain_samples([i])
        if isinstance(i, str):
            if match_prefix:
                self.samples.retain_samples_by_prefix([i])
            elif match_suffix:
                self.samples.retain_samples_by_suffix([i])
            else:
                self.samples.retain_samples_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            self.samples.retain_samples(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                self.samples.retain_samples_by_prefix(i)
            elif match_suffix:
                self.samples.retain_samples_by_suffix(i)
            else:
                self.samples.retain_samples_by_name(i)
        else:
            raise ValueError('i must be an int, str, list of int, or list of str.')

    # def insert_sites(self, sequence_str, i, marker_str=None):
    #     """Inserts a new sequence in the alignment matrix at the specified
    #     site position. This increases the total number of columns.

    #     Parameters
    #     ----------
    #     sequence_str : str or list of str
    #     i : int or list of int

    #     """
    #     if marker_str is None and self._marker_aln:
    #         assert ValueError('marker_str cannot be None if the alignment ' \
    #                           'has marker sequences')
    #     if (marker_str is not None) and (not self._marker_aln):
    #         assert ValueError('The alignment does not use marker sequences')
    #     self._sample_aln.insert_sites(sequence_str, i)
    #     if self._marker_aln:
    #         self._marker_aln.insert_sites(marker_str, i)

    # def append_site(self, sequence_str, marker_str=None):
    #     """Inserts a new sequence at after the last column of the
    #     alignment matrix. This increases the total number of columns by 1.

    #     Parameters
    #     ----------
    #     sequence_str : str

    #     """
    #     if marker_str is None and self._marker_aln:
    #         assert ValueError('marker_str cannot be None if the alignment ' \
    #                           'has marker sequences')
    #     if (marker_str is not None) and (not self._marker_aln):
    #         assert ValueError('The alignment does not use marker sequences')
    #     self._sample_aln.append_site(sequence_str)
    #     if self._marker_aln:
    #         self._marker_aln.append_site(marker_str)

    def remove_sites(self, i, description_encoder=None):
        """Removes sites based on the given index.
        If index is a number, only one site is removed.
        If the index is a list of numbers, the sequence found at each column
        number is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        if isinstance(i, int):
            i = [i]
        self.samples.remove_sites(i)
        if not (self.markers is None or self.markers.nsamples == 0):
            self.markers.remove_sites(i)
            assert self.samples.nsites == self.markers.nsites, \
                "Sample and marker nsites are not equal."
        # Update blocks if exists
        if self.blocklists:
            self.blocklists = [
                blockrs.remove_sites_from_blocks(blist, i)
                for seq, blist in zip(self.samples.sequences, self.blocklists)]
        if description_encoder:
            self.samples.set_descriptions(
                list(range(self.samples.nsamples)),
                [description_encoder(sid, blist)
                 for sid, blist in zip(self.samples.ids, self.blocklists)]
            )


    def retain_sites(self, i, description_encoder=None):
        """Keeps sites based on the given index.
        If index is a number, only one site is retained.
        If the index is a list of numbers, the characters at columns not
        found in the list is deleted.

        Parameters
        ----------
        i : int or list of int

        """
        if isinstance(i, int):
            i = [i]
        self.samples.retain_sites(i)
        if not (self.markers is None or self.markers.nsamples == 0):
            self.markers.retain_sites(i)
            assert self.samples.nsites == self.markers.nsites, \
                "Sample and marker nsites are not equal."
        # Update blocks if exists
        if self.blocklists:
            j = [pos for pos in range(self.samples.nsites) if pos not in i]
            self.blocklists = [
                blockrs.remove_sites_from_blocks(blist, j)
                for seq, blist in zip(self.samples.sequences, self.blocklists)]
        if description_encoder:
            self.samples.set_descriptions(
                list(range(self.samples.nsamples)),
                [description_encoder(sid, blist)
                 for sid, blist in zip(self.samples.ids, self.blocklists)]
            )


    def get_samples(self, i, match_prefix=False, match_suffix=False):
        """Returns a list of sequence strings containing only the samples
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        if isinstance(i, int):
            return self.samples.get_samples([i])
        elif isinstance(i, str):
            if match_prefix:
                return self.samples.get_samples_by_prefix([i])
            elif match_suffix:
                return self.samples.get_samples_by_suffixx([i])
            else:
                return self.samples.get_samples_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            return self.samples.get_samples(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            if match_prefix:
                return self.samples.get_samples_by_prefix(i)
            elif match_suffix:
                return self.samples.get_samples_by_suffixx(i)
            else:
                return self.samples.get_samples_by_name(i)
        else:
            raise ValueError('i must be an int, str, list of int, or list of str.')

    def get_markers(self, i):
        """Returns a list of sequence strings containing only the markers
        specified by the index.

        Parameters
        ----------
        i : int or list of int

        Returns
        -------
        list of str

        """
        if isinstance(i, int):
            return self.markers.get_samples([i])
        elif isinstance(i, str):
            return self.markers.get_samples_by_name([i])
        elif isinstance(i, list) and sum((isinstance(j, int) for j in i)):
            return self.markers.get_samples(i)
        elif isinstance(i, list) and sum((isinstance(j, str) for j in i)):
            return self.markers.get_samples_by_name(i)
        else:
            raise ValueError('i must be an int, str, list of int, or list of str.')

    # def get_sites(self, i):
    #     """Returns a new alignment containing only the sites specified
    #     by the index.

    #     Parameters
    #     ----------
    #     i : int or list of int

    #     Returns
    #     -------
    #     AlignmentMatrix

    #     """
    #     return self.__class__.subset(self, sites=i)

    @classmethod
    def from_fasta(cls, path, name, marker_kw=None):
        """Create an Alignment from a FASTA-formatted file.

        Parameters
        ----------
        path : str
            Path to FASTA file
        name : str
            name of alignment
        marker_kw : str, optional
            A sample is considered a marker if this keyword is present
            within the sequence ID

        Returns
        -------
        Alignment

        """
        d = fasta_file_to_lists(path, marker_kw=marker_kw)
        sample_aln = BaseAlignment(d['sample']['ids'],
                                   d['sample']['descriptions'],
                                   d['sample']['sequences'])
        marker_aln = BaseAlignment(d['marker']['ids'],
                                   d['marker']['descriptions'],
                                   d['marker']['sequences'])
        # Create alignments
        return cls(name, sample_aln, marker_aln)

    def to_fasta(self, path):
        """Saves the alignment as a FASTA-formatted text file.

        Parameters
        ----------
        path : str

        """
        with open(path, 'w') as writer:
            print(self.samples, file=writer)
            print(self.markers, file=writer)

    def set_blocklists(self, ref_seq, description_encoder=None):
        self.blocklists = [blockrs.pairwise_to_blocks(ref_seq, seq)
                           for seq in self.samples.sequences]
        if description_encoder:
            self.samples.set_descriptions(
                list(range(self.samples.nsamples)),
                [description_encoder(sid, blist)
                 for sid, blist in zip(self.samples.ids, self.blocklists)]
            )

    def parse_add_blocks(self, block_strings):
        if self.blocklists and len(block_strings) != self.samples.nsamples:
            raise ValueError('length of block string list not equal to the number of sequences')
        elif not self.blocklists:
            self.blocklists = [None for _ in range(self.samples.nsamples)]
        for i, block_str in enumerate(block_strings):
            # Parse block str into blocks
            blocks = []
            # Adds to blocklists
            self.blocklists[i] = blocks

    def __repr__(self):
        return '{}(nsamples={}, nsites={}, nmarkers={})'.format(
            self.__class__.__name__,
            self.nsamples,
            self.nsites,
            self.nmarkers
        )

    def __str__(self):
        return '\n'.join([str(self.samples), str(self.markers)])


# class CatBlock:
#     def __init__(self, name, start, stop):
#         self.name = name
#         self.start = start
#         self.stop = stop

#     def __repr__(self):
#         return 'CatBlock(name={}, start={}, stop={})'.format(
#             self.name, self.start, self.stop
#         )

#     def __str__(self):
#         return '{}={}:{}'.format(self.name, self.start, self.stop)

# class CatAlignment(Alignment):
#     """CatAlignment represents a superalignment of 2 or more
#     alignments concatenated together laterally.
#     """
#     def __init__(self, sequence_list, marker_list, catblocks, 
#                  block_lists_map=None, name=None,
#                  sample_to_uint_fn=None, uint_to_sample_fn=None,
#                  marker_to_uint_fn=None, uint_to_marker_fn=None):
#         """Creates a new CatAlignment.

#         Parameters
#         ----------
#         sequence_list : list of Sequence
#         marker_list : list of Marker
#         concat_list : list of CatBlock

#         """
#         super().__init__(sequence_list, marker_list, name=name,
#                          sample_to_uint_fn=sample_to_uint_fn,
#                          uint_to_sample_fn=uint_to_sample_fn,
#                          marker_to_uint_fn=marker_to_uint_fn,
#                          uint_to_marker_fn=uint_to_marker_fn)
#         self.catblocks = OrderedDict([(cb.id, self.catblocks[i])
#                                       for i, cb in enumerate(catblocks)])
#         self.block_lists_map = OrderedDict() if not block_lists_map else \
#                                block_lists_map

#     @classmethod
#     def concatenate(cls, aln_list, aln_ids=None, use_aln_names=True):
#         """Concatenates a list of Alignments into a single superalignment.

#         Parameters
#         ----------
#         aln_list : list of Alignment
#             [description]
#         aln_ids : list, optional
#             If specified, this list will be used as keys
#             to access individual alignments.
#         use_aln_names : bool, optional
#             When aln_ids is None, determines what values to use
#             as keys to access individual sequence alignments.
#             If True, alignment names are used . If False,
#             numbers from 0 corresponding to the position of the
#             alignment in the alignment list will be used.

#         Returns
#         -------
#         CatAlignment

#         """
#         start = 0
#         def coords(sid, val):
#             nonlocal start
#             start += val
#             return CatBlock(sid, start-val, start)
#         # Create a new concat alignment
#         new_aln = cls.__new__(cls)
#         new_aln.name = 'concat_' + '_'.join([str(aln.name) for aln in aln_list])
#         # Save alignment order
#         # Put block lists in a mapping
#         if aln_ids is not None:
#             new_aln.catblocks = OrderedDict([
#                 (i, coords(i, v.nsites)) for i, v in zip(aln_ids, aln_list)])
#             new_aln.block_lists_map = OrderedDict([
#                 (i, copy_block_lists(v.samples.block_lists))
#                 for i, v in zip(aln_ids, aln_list)])
#         elif use_aln_names:
#             new_aln.catblocks = OrderedDict([
#                 (v.name, coords(v.name, v.nsites)) for v in aln_list])
#             new_aln.block_lists_map = OrderedDict([
#                 (v.name, copy_block_lists(v.samples.block_lists))
#                 for v in aln_list])
#         else:
#             new_aln.catblocks = OrderedDict([
#                 (i, coords(i, v.nsites)) for i, v in enumerate(aln_list)])
#             new_aln.block_lists_map = OrderedDict([
#                 (i, copy_block_lists(v.samples.block_lists))
#                 for i, v in enumerate(aln_list)])
#         # Create new concatenated block lists
#         total_sites = sum((aln.nsites for aln in aln_list))
#         concat_block_lists = [[Block(0, total_sites)]
#                               for i in range(aln_list[0].nsamples)]
#         # Create new sample alignment from matrix
#         new_aln._sample_aln = SampleAlignment.from_uint_matrix(
#             np.concatenate([aln.sample_matrix for aln in aln_list], axis=1),
#             aln_list[0].samples.ids,
#             [','.join([str(aln.name) for aln in aln_list])] * \
#                 len(aln_list[0].samples.descriptions),  # replaces desc
#             concat_block_lists,  # empties block list
#             to_uint_fn=aln_list[0].samples.custom_to_uint_fn,
#             from_uint_fn=aln_list[0].samples.custom_from_uint_fn,
#             to_block_fn=aln_list[0].samples.custom_to_block_fn,
#             from_block_fn=aln_list[0].samples.custom_from_block_fn,
#         )
#         # Create new marker alignment from matrix
#         new_aln._marker_aln = MarkerAlignment.from_uint_matrix(
#             np.concatenate([aln.marker_matrix for aln in aln_list], axis=1),
#             aln_list[0].markers.ids,
#             aln_list[1].markers.descriptions,
#             to_uint_fn=aln_list[0].samples.custom_to_uint_fn,
#             from_uint_fn=aln_list[0].samples.custom_from_uint_fn,
#         )
#         return new_aln

#     @classmethod
#     def from_fasta(cls, path, name=None, marker_kw=None,
#                    block_lists_map=None,
#                    sample_to_uint_fn=None, uint_to_sample_fn=None,
#                    marker_to_uint_fn=None, uint_to_marker_fn=None):
#         """Create a CatAlignment from a FASTA-formatted file.

#         Parameters
#         ----------
#         path : str
#             Path to FASTA file
#         marker_kw : str, optional
#             A sample is considered a marker if this keyword is present
#             within the sequence ID
#         sample_to_uint_fn : function, optional
#         uint_to_sample_fn : function, optional
#         marker_to_uint_fn : function, optional
#         uint_to_marker_fn : function, optional

#         Raises
#         ------
#         TypeError

#         Returns
#         -------
#         Alignment

#         """
#         sequence_list = []
#         marker_list = []
#         for item in fasta_file_to_list(path, marker_kw=marker_kw):
#             if isinstance(item, Sequence):
#                 sequence_list.append(item)
#             elif isinstance(item, Marker):
#                 marker_list.append(item)
#             else:
#                 raise TypeError('expected Sequence or Marker object')
#         catblocks = string_to_catblocks(sequence_list[0].description.rstrip())
#         return cls(sequence_list, marker_list, catblocks, name=name,
#                    block_lists_map=block_lists_map,
#                    sample_to_uint_fn=sample_to_uint_fn,
#                    uint_to_sample_fn=uint_to_sample_fn,
#                    marker_to_uint_fn=marker_to_uint_fn,
#                    uint_to_marker_fn=uint_to_marker_fn)

#     def to_fasta(self, path, catblocks_path=None, block_lists_path=None):
#         """Saves the concatenated alignment as a FASTA-formatted text file.
#         Catblocks and block lists can also be simultaneously saves as
#         tab-delimitted files.

#         Parameters
#         ----------
#         path : str
#         catblocks_path : str, optional
#         block_lists_path : str, optional

#         """
#         # Overwrite descriptions with catblocks
#         self.samples.descriptions = [catblocks_to_string(self.catblocks)
#                                      for i in self.nsamples]
#         if catblocks_path:
#             with open(catblocks_path, 'w') as cb_writer:
#                 print('{}\t{}\t{}'.format('name', 'start', 'stop'), 
#                       file=cb_writer)
#                 for cb in self.catblocks.values():
#                     print('{}\t{}\t{}'.format(cb.name, cb.start, cb.stop), 
#                           file=cb_writer)
#         if block_lists_path:
#             with open(block_lists_path, 'w') as b_writer:
#                 print('{}\t{}\t{}'.format('name', 'start', 'stop'), 
#                       file=b_writer)
#                 for name, blist in self.block_lists_map.items():
#                     for b in blist:
#                         print('{}\t{}\t{}'.format(name, b.start, b.stop),
#                               file=b_writer)
#         with open(path, 'w') as writer:
#             print(self, file=writer)

#     def get_alignment(self, name):
#         """Return the subalignment by name or key.

#         Parameters
#         ----------
#         name : int or str
#             Name assigned to the alignment during concatenation.

#         Raises
#         ------
#         IndexError
#             Returns an IndexError when the given name does not
#             match any alignment key.

#         Returns
#         -------
#         Alignment

#         """

#         if name not in self.catblocks:
#             raise IndexError("name not found")
#         _, start, stop = self.catblocks[name]
#         aln = Alignment.subset(self, sites=range(start, stop))
#         aln.name = 'subaln_{}'.format(name)
#         aln.samples.block_lists = copy_block_lists(self.block_lists_map[name])
#         return aln

#     def splitg(self):
#         """Splits the concatenated superalignment into individual subalignments
#         and returns a generator.

#         Yields
#         ------
#         Alignment

#         """

#         return (Alignment.subset(self, sites=range(cb.start, cb.stop))
#                 for cb in self.catblocks)

#     def split(self):
#         """Splits the concatenated superalignment as a list of its
#         individual subalignments.

#         Returns
#         -------
#         list of Alignment

#         """
#         return list(self.splitg())


# def fasta_file_to_list(path, marker_kw=None):
#     """Reads a FASTA formatted text file to a list.

#     Parameters
#     ----------
#     path : str

#     Returns
#     -------
#     list of tuple

#     """
#     name = ''
#     description = ''
#     _seq = ''
#     seq_list = []
#     with open(path, 'r') as f:  # pylint: disable=invalid-name
#         for line in f.readlines():
#             line = line.rstrip()
#             if line.startswith('>'):
#                 # Store sequence if _seq has contents
#                 if _seq:
#                     if marker_kw:
#                         if marker_kw in name:
#                             seq = Marker(name, description, _seq)
#                         else:
#                             seq = Sequence(name, description, _seq)
#                     else:
#                         seq = Sequence(name, description, _seq)
#                     seq_list.append(seq)
#                     _seq = ''
#                 # Split id and description
#                 try:
#                     name, description = line[1:].split(' ', 1)
#                 except ValueError:
#                     name = line[1:]
#                     description = ''
#             else:
#                 _seq += line
#         if _seq:
#             if marker_kw:
#                 if marker_kw in name:
#                     seq = Marker(name, description, _seq)
#                 else:
#                     seq = Sequence(name, description, _seq)
#             else:
#                 seq = Sequence(name, description, _seq)
#             seq_list.append(seq)
#     return seq_list

def fasta_file_to_lists(path, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str
    name : str

    Returns
    -------
    Alignment

    """
    _id = ''
    _description = ''
    _seq = ''

    sample_ids = []
    sample_descs = []
    sample_seqs = []

    marker_ids = []
    marker_descs = []
    marker_seqs = []

    with open(path, 'r') as f:  # pylint: disable=invalid-name
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                # Store sequence if _seq has contents
                if _seq:
                    if marker_kw and (marker_kw in _id):
                        marker_ids.append(_id)
                        marker_descs.append(_description)
                        marker_seqs.append(_seq)
                    else:
                        sample_ids.append(_id)
                        sample_descs.append(_description)
                        sample_seqs.append(_seq)
                    _seq = ''
                # Split id and description
                try:
                    _id, _description = line[1:].split(' ', 1)
                except ValueError:
                    _id, _description = line[1:], ''
            else:
                _seq += line
        if _seq:
            if marker_kw and (marker_kw in _id):
                marker_ids.append(_id)
                marker_descs.append(_description)
                marker_seqs.append(_seq)
            else:
                sample_ids.append(_id)
                sample_descs.append(_description)
                sample_seqs.append(_seq)
    return {
        'sample': {
            'ids': sample_ids,
            'descriptions': sample_descs,
            'sequences': sample_seqs,
        },
        'marker': {
            'ids': marker_ids,
            'descriptions': marker_descs,
            'sequences': marker_seqs,
        }
    }

def fasta_file_to_alignment(path, name, marker_kw=None):
    """Reads a FASTA formatted text file to a list.

    Parameters
    ----------
    path : str
    name : str

    Returns
    -------
    Alignment

    """
    d = fasta_file_to_lists(path, marker_kw=marker_kw)
    sample_aln = BaseAlignment(d['sample']['ids'],
                               d['sample']['descriptions'],
                               d['sample']['sequences'])
    marker_aln = BaseAlignment(d['marker']['ids'],
                               d['marker']['descriptions'],
                               d['marker']['sequences'])
    # Create alignments
    return Alignment(name, sample_aln, marker_aln)


# def copy_block_lists(block_lists):
#     return [[Block(b.start, b.stop) for b in blist] for blist in block_lists]

# def catblocks_to_string(catblock_list):
#     return ';'.join([str(cb) for cb in catblock_list])

# def string_to_catblocks(string, int_names=False):
#     catblocks_raw = re.findall(r'(\S+?)\=(\d+?)\:(\d+?)', string)
#     if int_names:
#         return [CatBlock(int(i[0]), int(i[1]), int(i[2]))
#                 for i in catblocks_raw]
#     return [CatBlock(i[0], int(i[1]), int(i[2])) for i in catblocks_raw]
