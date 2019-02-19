from copy import deepcopy

from libalignmentrs.record import Record

# Shared functions

# Insert Methods
# ==========================================================================


# Insert/append many samples
# ------------------------------
# TODO: Use record for insert/append many
def insert_items(baln, row, records, copy=False):
    raise NotImplementedError()

def append_items(baln, records, copy=False):
    raise NotImplementedError()

def insert_items_from_lists(baln, row, ids, descriptions, sequences):
    """Inserts sample entries at the specified position.

    Parameters
    ----------
    row : int
        Row position to insert the samples.
    ids : list of str
        List of sample identifiers to insert. The order should correspond
        to the order of ther other lists.
    descriptions : list of str
        List of sample descriptions to insert. The order should correspond
        to the order of ther other lists.
    sequences : list of str
        List of sample sequneces to insert. The order should correspond
        to the order of ther other lists.

    """
    if not(isinstance(ids, list) and
            sum((isinstance(j, str) for j in ids))):
        raise TypeError('ids must be a list of str.')
    if not(isinstance(descriptions, list) and
            sum((isinstance(j, str) for j in descriptions))):
        raise TypeError('descriptions must be a list of str.')
    if not(isinstance(sequences, list) and
            sum((isinstance(j, str) for j in sequences))):
        raise TypeError('sequences must be a list of str.')
    baln.insert_rows(row, ids, descriptions, sequences)


def append_items_from_lists(baln, ids, descriptions, sequences):
    """Appends sample entries after the last sample in the alignment.

    Parameters
    ----------
    rows : list of str
        List of sample identifiers to insert. The order should correspond
        to the order of ther other lists.
    descriptions : list of str
        List of sample descriptions to insert. The order should correspond
        to the order of ther other lists.
    sequences : list of str
        List of sample sequneces to insert. The order should correspond
        to the order of ther other lists.

    """
    if not(isinstance(ids, list) and
            sum((isinstance(j, str) for j in ids))):
        raise TypeError('ids must be a list of str.')
    if not(isinstance(descriptions, list) and
            sum((isinstance(j, str) for j in descriptions))):
        raise TypeError('descriptions must be a list of str.')
    if not(isinstance(sequences, list) and
            sum((isinstance(j, str) for j in sequences))):
        raise TypeError('sequences must be a list of str.')
    baln.append_rows(ids, descriptions, sequences)


# Setter methods
# ==========================================================================
def replace_items(baln, rows, sequences):
    """Replaces sequences of the specified items.

    Parameters
    ----------
    rows : int, str, list of int, or list of str
        An int/str/list specifying the items to be replaced.
        If a list, length must match the number of given sequences.
    sequences : str
        A single sequence or a list of sequences.
        If `rows` is an int/str, `sequences` must be a str.
        If `rows` is a list, `sequences` must list and the number of sequences
        must match the length of `rows`.

    """
    # Calls specific replace_sequence setter depending on the
    # type if rows
    if isinstance(rows, int) and isinstance(sequences, str):
        baln.replace_sequences([rows], [sequences])
    elif isinstance(rows, str) and isinstance(sequences, str):
        rows = baln.row_names_to_ids([rows])
        # TODO: Change to error
        assert len(rows) == 1, \
            'Length of matched rows and number of provided ' \
            'sequences do not match.'
        baln.replace_sequences([rows], [sequences])
    elif isinstance(rows, list) and sum((isinstance(j, int) for j in rows)):
        baln.replace_sequences(rows, sequences)
    elif isinstance(rows, list) and sum((isinstance(j, str) for j in rows)):
        rows = baln.row_names_to_ids(rows)
        # TODO: Change to error
        assert len(rows) == len(sequences), \
            'Length of matched rows and number of provided ' \
            'sequences do not match.'
        baln.replace_sequences(rows, sequences)
    else:
        raise TypeError('rows must be an int, str, list of int, or list of str.')


def reorder_items(baln, rows):
    """Reorders items based on a list of identifiers or indices.

    Parameters
    ----------
    rows : list of int or list of str
        Order of the list specifies the new ordering of the items.

    """
    # Check type of rows, and convert if necessary
    if isinstance(rows, list) and sum((isinstance(j, int) for j in rows)):
        pass
    elif isinstance(rows, list) and sum((isinstance(j, str) for j in rows)):
        rows = baln.row_names_to_ids(rows)
    else:
        raise TypeError('rows must be a list of int or list of str.')
    baln.reorder_rows(rows)


# deleters
# ------------------------------
def remove_items(baln, rows, match_prefix=False, match_suffix=False):
    """Removes items based a list of identifiers or indices.
    This is the opposite of the `retain_items` method.

    Parameters
    ----------
    rows : int, str, list of int, or list of str
        An int/str/list specifying the items to be removed.
    match_prefix : bool, optional
        Whether to interpret `rows` as a prefix to match against
        the list of marker names. (default is False)
    match_suffix : bool, optional
        Whether to interpret `rows` as a suffix to match against
        the list of marker names. This parameter is considered
        only if match_prefix is False. (default is False)

    """
    if isinstance(rows, int):
        baln.remove_rows([rows])
    elif isinstance(rows, str):
        if match_prefix:
            baln.remove_rows_by_prefix([rows])
        elif match_suffix:
            baln.remove_rows_by_suffix([rows])
        else:
            baln.remove_rows_by_name([rows])
    elif isinstance(rows, list) and sum((isinstance(j, int) for j in rows)):
        baln.remove_rows(rows)
    elif isinstance(rows, list) and sum((isinstance(j, str) for j in rows)):
        if match_prefix:
            baln.remove_rows_by_prefix(rows)
        elif match_suffix:
            baln.remove_rows_by_suffix(rows)
        else:
            baln.remove_rows_by_name(rows)
    else:
        raise TypeError('rows must be an int, str, list of int, or list of str.')


def retain_items(baln, rows, match_prefix=False, match_suffix=False):
    """Removes items based a list of identifiers or indices.
    This is the opposite of the `remove_items` method.

    Parameters
    ----------
    rows : int, str, list of int, or list of str
        An int/str/list specifying the items to be retained.
    match_prefix : bool, optional
        Whether to interpret `rows` as a prefix to match against
        the list of marker names. (default is False)
    match_suffix : bool, optional
        Whether to interpret `rows` as a suffix to match against
        the list of marker names. This parameter is considered
        only if match_prefix is False. (default is False)

    """
    if isinstance(rows, int):
        baln.retain_rows([rows])
    elif isinstance(rows, str):
        if match_prefix:
            baln.retain_rows_by_prefix([rows])
        elif match_suffix:
            baln.retain_rows_by_suffix([rows])
        else:
            baln.retain_rows_by_name([rows])
    elif isinstance(rows, list) and sum((isinstance(j, int) for j in rows)):
        baln.retain_rows(rows)
    elif isinstance(rows, list) and sum((isinstance(j, str) for j in rows)):
        if match_prefix:
            baln.retain_rows_by_prefix(rows)
        elif match_suffix:
            baln.retain_rows_by_suffix(rows)
        else:
            baln.retain_rows_by_name(rows)
    else:
        raise TypeError('rows must be an int, str, list of int, or list of str.')


# Iterators
# ==========================================================================
def iter_aln_sites(baln, start=0, stop=None, size=1):
    """Iterates column-wise over the base alignment.

    Parameters
    ----------
    start : int, optional
        Starting position. (default is 0)
    stop : [type], optional
        Stopping column position.
        If None, the iterator will continue until the end of the alignment.
        (default is None)
    size : int, optional
        Number of characters to yield at each iteration.
        For single characters, `size` = 1.
        For codons, `size` = 3. (default is 1)

    Raises
    ------
    ValueError
        If the alignment cannot be cleanly cut up into the specified
        chunk size (`ncols` not divisible be `size`),
        a ValueError is raised.

    Yields
    ------
    list of str
        List of sequences representing a site or chunk of the alignment.

    """
    if stop is None:
        stop = baln.ncols
    if (stop - start) % size != 0:
        raise ValueError('Alignment cannot be completely divided into '
                        'chucks of size {}'.format(size))
    for i in range(start, stop, size):
        yield [s[i:i+size] for s in baln.sequences]

def iter_aln(baln):
    """Iterates over alns in the alignment.
    Excludes alns.

    Yields
    ------
    str

    """
    for sequence in baln.sequences:
        yield sequence

def iter_aln_records(baln):
    """Iterates over alns in the alignment, returning a Record object.
    Excludes alns.

    Yields
    ------
    Record

    """
    for vals in zip(baln.ids,
                    baln.descriptions,
                    baln.sequences):
        yield Record(*vals)



def subset(aln, cols, **kwargs):
        """Returns a subset of a given alignment based on the
        specified set of samples, markers and sites.

        Parameters
        ----------
        aln : Alignment
        cols : int, list of int, or None
            An int/list specifying the sites to be included.
            If None, all sites will be included in the subset.
        **kwargs
            An int/str/list specifying the samples to be included.
            If None, all samples will be included in the subset.

        Raises
        ------
        TypeError
            Given parameter has the wrong parameter type.
        ValueError
            marker_ids is specified by aln.markers is empty.

        Returns
        -------
        Alignment
            New alignment object containing the subset of sample and
            markers rows, and site columns.
            This subset is a deep copy of the original alignment and
            will not be affect by changes made in the original.

        """
        # Checks the value of sample_ids and converts if necessary.
        new_alns = []
        ro = 0
        # Check if keywords in kwargs are all members
        for key in kwargs.keys():
            if key not in aln.__class__.members:
                raise TypeError(
                    'Member not found: {} is an invalid keyword.'.format(key))

        # Checks the value of cols and converts if necessary.
        if cols is None:
            cols = list(range(0, aln.ncols))
        elif isinstance(cols, int):
            cols = [cols]
        elif (isinstance(cols, list) and
              sum((isinstance(j, int) for j in cols))):
            pass
        else:
            raise TypeError('cols must be an int, or list of int.')

        for member in aln.__class__.members:
            if member in kwargs.keys():
                rows = kwargs[member]
                # Checks the value of rows and converts if necessary.
                if rows is None:
                    rows = list(range(0, self.nrows))
                elif isinstance(rows, int):
                    rows = [rows]
                elif isinstance(rows, str):
                    rows = aln.__getattribute__(member).row_names_to_ids([rows])
                elif (isinstance(rows, list) and
                    sum((isinstance(j, int) for j in rows))):
                    pass
                elif (isinstance(rows, list) and
                    sum((isinstance(j, str) for j in rows))):
                    rows = aln.__getattribute__(member).row_names_to_ids(rows)
                else:
                    raise TypeError('rows must be an int, str, list of int, '
                                    'or list of str.')
                rows = [row-ro for row in rows if row-ro > 0]
                if rows:
                    new_aln = aln.__getattribute__(member).subset(rows, cols)
                else:
                    new_aln = None
            else:
                new_aln = aln.__getattribute__(member).copy()
            new_alns.append(new_aln)
            ro += aln.__getattribute__(member).nrows
        
        return aln.__class__(
            aln.name, *new_aln,
            linspace=aln.__getattribute__('_linspace').extract(cols),
            metadata=deepcopy(aln.metadata))
