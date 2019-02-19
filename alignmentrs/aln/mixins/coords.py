import re
from libalignmentrs.position import Block, BlockSpace


class CoordsMixin:
    _description_re = re.compile(r'\s*coords\=\{.+\}\s*')

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

        new_descriptions = [regex.sub(' ', desc) for desc in self.descriptions]
        self.replace_descriptions(
            list(range(len(self.descriptions))), new_descriptions)

    @staticmethod
    def _default_description_encoder(current_description: str, blocks: list):
        current_description = \
            CoordsMixin._description_re.sub(' ', current_description)
        blocks_str = \
            'coords={' + \
            ';'.join(['{}:{}'.format(b.start, b.stop) for b in blocks]) + \
            '}'
        return ' '.join([current_description, blocks_str])
