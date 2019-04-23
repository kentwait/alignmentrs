

class MockData:
    def __init__(self):
        self.data = [
            'ATGCAT',
            'ATGGGT',
            'ATGAAT',
        ]

    def get_row(self, i):
        return self.sequences[i]