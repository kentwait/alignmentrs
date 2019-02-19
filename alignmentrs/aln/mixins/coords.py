from libalignmentrs.position import BlockSpace


class CoordsMixin:
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