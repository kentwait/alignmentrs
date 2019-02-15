import numpy as np


def mark_sites_with_chars(aln, target_list, size=1,
                          ignore_case=True, inverse=False,
                          copy=False, _sep='|'):
    """Adds markers to the alignment indicating which sites
    matche/contain the target character or string.

    Adds one marker for each item in `target_list`.

    Parameters
    ----------
    aln : Alignment
        Alignment to screen.
    target_list : list of str
        List of target characters (ie. 'N' for ambiguous characters
        or '-' for gaps). Each target generates a marker track.
        If a target is itself a list, then each site will be tested
        if it is one of the members of the list.
    size : int, optional
        Sie of the size in terms of number of alignment columns.
        For single characters such as nucleotides, size = 1.
        For codons, size = 3. (default is 1,)
    ignore_case : bool, optional
        Whether to consider upper and lower case letters to be the same.
        (the default is True)
    inverse : bool, optional
        When `inverse` is True, matching sites are marked with 0 while
        non-matching sites are marked with 1. (the default is False,
        which marks matching site 1, otherwise 0)
    copy : bool, optional
        Returns a new copy instead of adding new markers inplace.
        (default is False, operation is done inplace)

    Raises
    ------
    ValueError
        When the number of columns in the alignment is not divisible by
        the specified size, a ValueError is raised.

    Returns
    -------
    Alignment or None
        If copy is True, returns a new alignment, otherwise no
        value is returned (None).

    """
    aln = aln.__class__(
        aln.name, aln.samples.copy(), aln.markers.copy()) if copy else \
        aln
    if aln.nsites % size != 0:
        raise  ValueError('Alignment cannot be completely divided into '
                          'chucks of size {}'.format(size))

    position_list = []
    changer = lambda x: x.upper() if ignore_case else x
    t_c, f_c = ('0', '1') if inverse else ('1', '0')
    for target in target_list:
        # Create an initial filter array of 1
        filter_array = np.ones(int(aln.nsites/size))

        # Determine sites with char in within the site
        if isinstance(target, list):
            position_list = [
                i
                # Loop over sample sites by size steps,
                # Sites is a list of size-char strings
                for i, sites in enumerate(aln.iter_sample_sites(size=size))
                # Loop over each unique variant of strings
                for variant in set(sites)
                # If target is found, include the current position i
                if changer(variant) in [changer(t) for t in target]
            ]
            target_name = _sep.join(target)
        else:
            position_list = [
                i
                # Loop over sample sites by size steps,
                # Sites is a list of size-char strings
                for i, sites in enumerate(aln.iter_sample_sites(size=size))
                # Loop over each unique variant of strings
                for variant in set(sites)
                # If target is found, include the current position i
                if changer(target) in changer(variant)
            ]
            target_name = target
        filter_array[position_list] = 0

        # Add new marker
        aln.markers.append_rows(
            ['{}_marker'.format(target_name)],
            ['notes="{} if site has "{}", else {}"'.format(
                t_c*size, target, f_c*size)],
            [''.join([t_c*size if i else f_c*size for i in filter_array])]
        )
    if copy:
        return aln


def drop_sites_using_binary_markers(aln, marker_ids, inverse=False,
                                    match_prefix=False, match_suffix=False,
                                    copy=False):
    """Removes sites that failed to pass all of the markers in a given
    list of markers.

    This function assumes that the selected markers are binary-valued
    markers whose values are "0" and "1". A column passed a marker if the
    column contains a "1", otherwise it failed.

    Parameters
    ----------
    aln : Alignment
        Alignment containing the markers to remove sites.
    marker_ids : int, str, list of int, or list of str
        int, str or list specifying the markers to be used.
    inverse : bool, optional
        When inverse if True, columns that passed all the selected markers
        (all 1's) will be removed. (default is False, columns that passed
        the selected markers will be kept)
    match_prefix : bool, optional
        Whether to interpret `i` as a prefix to match against
        the list of markers names. (default is False)
    match_suffix : bool, optional
        Whether to interpret `i` as a suffix to match against
        the list of markers names. This parameter is considered
        only if match_prefix is False. (default is False)
    copy : bool, optional
        Returns a new copy instead of performing dropping inplace.
        (default is False, operation is done inplace)

    Returns
    -------
    Alignment or None
        If copy is True, returns a new alignment, otherwise no
        value is returned (None).

    """
    aln = aln.copy() if copy else aln
    # Get marker alignments and turn into a numpy array
    marker_matrix = np.array(
        [list(map(int, m))
         for m in aln.get_markers(marker_ids,
                                  match_prefix=match_prefix,
                                  match_suffix=match_suffix)
         .sequences])
    # Sum the values down each column
    # Columns whose sum is less than the number of rows have failed
    # one or more filters
    summed = np.sum(marker_matrix, axis=0)
    remove_list = np.where(summed == len(marker_matrix))[0] if inverse else \
                  np.where(summed < len(marker_matrix))[0]

    # Edit alignment inplace
    aln.remove_sites(remove_list)

    if copy:
        return aln

def aln_to_sample_matrix(aln, size=1):
    """Convert an alignment's sample sequences into a numpy matrix.

    Parameters
    ----------
    size : int, optional
        Defines the number of characters is in each cell.
        For example, for single characters, set `size` = 1,
        while for a codon-based matrix, set `size` = 3.

    Returns
    -------
    np.array

    """
    return np.array([list(s) for s in aln.iter_sample_sites(size=size)]).T

def aln_to_marker_matrix(aln, size=1):
    """Converts an alignment's marker sequences into a numpy matrix.

    Parameters
    ----------
    size : int, optional
        Defines the number of characters is in each cell.
        For example, for single characters, set `size` = 1,
        while for a codon-based matrix, set `size` = 3.

    Returns
    -------
    numpy.array
        The multiple sequence alignment is converted into a numpy matrix
        with a shape corresponding to the number of samples and sites,
        respectively.

    """
    return np.array([list(s) for s in aln.iter_marker_sites(size=size)]).T
