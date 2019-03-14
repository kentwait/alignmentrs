

class MockData:
    def __init__(self):
        self.sequences = [
            'ATGCAT',
            'ATGGGT',
            'ATGAAT',
        ]

    def get_row(self, i):
        return self.sequences[i]