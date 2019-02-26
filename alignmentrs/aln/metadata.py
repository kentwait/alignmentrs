
class MetadataRedirect:
    def __init__(self, instance):
        self._instance = instance
        self._axis = 0

    @property
    def rows(self):
        """pandas.core.frame.DataFrame: Returns the associated row
        metadata as a pandas DataFrame."""
        return self._instance._row_metadata

    @property
    def cols(self):
        """pandas.core.frame.DataFrame: Returns the associated column
        metadata as a pandas DataFrame."""
        return self._instance._column_metadata

    @property
    def comments(self):
        """dict: Returns the comments for this alignment."""
        return self._instance._comments



