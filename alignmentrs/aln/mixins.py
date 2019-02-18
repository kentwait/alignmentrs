import itertools


class PropsMixin:
    @property
    def nrows(self):
        """int: Returns the number of rows in the alignment."""
        return sum(
            (self.__getattribute__(m).nrows for m in self.__class__.members))

    @property
    def ncols(self):
        """int: Returns the number of sites in the alignment."""
        return self.samples.ncols

    @property
    def ids(self):
        """list of str: Returns the list of identifiers."""
        generator = (
            self.__getattribute__(m).ids for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def descriptions(self):
        """list of str: Returns the list of descriptions."""
        generator = (
            self.__getattribute__(m).descriptions 
            for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def sequences(self):
        """list of str: Returns the list of sequences."""
        generator = (
            self.__getattribute__(m).sequences 
            for m in self.__class__.members)
        return [item for item in itertools.chain(*generator)]

    @property
    def coordinates(self):
        """list of int: Returns the list of coordinates
        associated to each column in the alignment."""
        return self._linspace.to_arrays()[0]


class SamplePropsMixin:
    @property
    def nsamples(self):
        """int: Returns the number of samples in the alignment."""
        if not self.markers:
            return []
        if not self.samples:
            return 0
        return self.samples.nrows

    @property
    def sample_ids(self):
        """list of str: Returns the list of sample identifiers."""
        if not self.samples:
            return []
        return self.samples.ids

    @property
    def sample_descriptions(self):
        """list of str: Returns the list of sample descriptions."""
        if not self.samples:
            return []
        return self.samples.descriptions

    @property
    def sample_sequences(self):
        """list of str: Returns the list of sample sequences."""
        if not self.samples:
            return []
        return self.samples.sequences


class MarkerPropsMixin:
    @property
    def nmarkers(self):
        """int: Returns the number of markers in the alignment."""
        if not self.markers:
            return 0
        return self.markers.nrows

    @property
    def marker_ids(self):
        """list of str: Returns the list of markers identifiers."""
        if not self.markers:
            return []
        return self.markers.ids

    @property
    def marker_descriptions(self):
        """list of str: Returns the list of marker descriptions."""
        if not self.markers:
            return []
        return self.markers.descriptions

    @property
    def marker_sequences(self):
        """list of str: Returns the list of marker sequences."""
        if not self.markers:
            return []
        return self.markers.sequences


class SampleAlnMixin:
    pass


class MarkerAlnMixin:
    pass
