import warnings


__all__ = ['NoNameWarning', 'DuplicateNameWarning']


class NoNameWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass

class DuplicateNameWarning(UserWarning):
    """Warning for mismatched/incompatible alignments.
    """
    pass