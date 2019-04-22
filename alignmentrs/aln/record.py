import pandas

class Record:
    def __init__(self, record_id, sequence, metadata):
        self.id = record_id
        self.sequence = self._make_sequence(sequence)
        self.metadata = self._make_metadata(metadata)

    def _make_sequence(self, sequence=None):
        if sequence is None:
            return ''
        elif isinstance(sequence, str):
            return sequence
        elif (isinstance(sequence, list) or isinstance(sequence, tuple)):
            try:
                string = ''.join(sequence)
            except TypeError:
                raise TypeError('Cannot construct `sequence` using the given `sequence` input: {}'.format(sequence))
        # Unsupported data type
        raise TypeError('`sequence` must be a str, list or tuple, instead got: {}'.format(type(sequence)))

    def _make_metadata(self, metadata=None):
        if metadata is None:
            return pandas.DataFrame(None)
        elif isinstance(metadata, pandas.DataFrame):
            return metadata
        # Unsupported data type
        raise TypeError('`metadata` must be a pandas.DataFrame, instead got: {}'.format(type(metadata)))
