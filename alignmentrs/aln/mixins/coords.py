import re
from libalignmentrs.position import Block, BlockSpace
from libalignmentrs.position import simple_block_str_to_linspace


class CoordsMixin:
    _description_re = re.compile(r'\s*coords\=\{(.+)\}\s*')

    @property
    def coordinates(self):
        """list of int: Returns the list of coordinates
        associated to each column in the alignment."""
        return self._linspace.to_arrays()[0]
    
    def reset_coordinates(self, start=0, stop=None, state=1):
        """Resets the coordinates of the alignment columns.  

        Parameters
        ----------
        start : int, optional
            Starting coordinate value (default is 0)
        stop : int, optional
            End coordinate value.
            This value is not part of linear space representing the alignment.
            (default is None, ending value becomes the sum of start and
            the number of sites)
        state : int, optional
            State annotation. (default is 1)

        """
        start = start if start is not None else 0
        stop = stop if stop is not None else start + self.samples.ncols
        state = state if state is not None else 1
        self._linspace: BlockSpace = BlockSpace(start, stop, state)

    def reset_coords_from_description(self, description_decoder: callable=None,
                                      read_from: int=None):
        if description_decoder is None:
            description_decoder = self._default_description_decoder
        if read_from is not None:
            if isinstance(read_from, int):
                if read_from >= self.nrows:
                    raise IndexError('read_from value is invalid.')
                read_from_ids = [read_from]
            else:
                raise ValueError('read_from must be an int.')
        else:
            read_from_ids = range(self.nrows)
        
        linspace = None
        for i in read_from_ids:
            linspace = description_decoder(self.descriptions[i])
            # Break once a linspace was sucessfully read
            if linspace is not None:
                break

        if linspace is None:
            raise ValueError(
                'Coordinates could not be parsed from the description')
        # Make a setter to check type and if length == ncols
        self._linspace = linspace

    @staticmethod
    def _default_description_decoder(description):
        match = CoordsMixin._description_re.search(description)
        try:
            linspace = simple_block_str_to_linspace(match.group(1))
        except AttributeError:
            linspace = None
        return linspace

    def write_coords_to_description(self, description_encoder: callable=None,
                                    write_to: int=None):
        if description_encoder is None:
            description_encoder = self._default_description_encoder
        if write_to is not None:
            if isinstance(write_to, int):
                if write_to >= self.nrows:
                    raise IndexError('write_to value is invalid.')
                write_to_ids = [write_to]
            else:
                raise ValueError('write_to must be an int.')
        else:
            write_to_ids = range(self.nrows)
        
        blocks = self._linspace.to_blocks()
        for i in write_to_ids:
            new_description = description_encoder(self.descriptions[i], blocks)
            self.replace_descriptions([i], [new_description])

    def remove_coords_from_description(self, coords_pattern=None):
        new_descriptions = []
        if coords_pattern is None:
            regex = CoordsMixin._description_re
        else:
            regex = re.compile(coords_pattern)

        new_descriptions = [regex.sub(' ', desc).strip()
                            for desc in self.descriptions]
        self.replace_descriptions(
            list(range(len(self.descriptions))), new_descriptions)

    @staticmethod
    def _default_description_encoder(current_description: str, blocks: list):
        parts = []
        blocks_str = \
            'coords={' + \
            ';'.join(['{}:{}'.format(b.start, b.stop) for b in blocks]) + \
            '}'
        parts.append(blocks_str)

        if current_description:
            current_description = \
                CoordsMixin._description_re.sub(' ', current_description)
            parts.append(current_description)
            
        return ' '.join(parts)
