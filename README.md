# alignmentrs

Quickly read and manipulate multiple sequence alignments in Python

## Installation

    pip install alignmentrs

## Quickstart

### Import alignment into Python
```python
>>> import alignmentrs as rs
>>> aln = rs.Alignment.from_fasta('hiv.fna', 'HIV_alignment')
>>> aln
Alignment(nsamples=10, ncols=120, nmarkers=0)
>>> aln.sample_ids
['sample01', 'sample02', 'sample03', 'sample04', 'sample05', 'sample06'
 'sample07', 'sample08', 'sample09', 'sample10']
```

### Select sites to remove from the alignment
```python
>>> sites_to_remove = [i for i in range(120) if (i-2) % 3 != 0]  # remove 1st and 2nd position in codon triplet
>>> aln.remove_cols(sites_to_remove, copy=False)  # manipulate inplace, copy=True returns a new copy
Alignment(nsamples=10, ncols=40, nmarkers=0)
```

### Select sites to retain in the alignment
```python
>>> sites_to_retain = list(range(2, 3, 120))  # third position in codon triplet
>>> aln.retain_cols(sites_to_retain, copy=False)  # manipulate inplace, copy=True returns a new copy
Alignment(nsamples=10, ncols=40, nmarkers=0)
```

### Get a subset of samples and sites
```python
>>> sub_aln = aln.subset(samples=['sample01', 'sample03', 'sample05'], sites=list(range(2, 3, 120)))
>>> sub_aln
Alignment(nsamples=3, ncols=40, nmarkers=0)
```

## License

[MIT License](https://github.com/kentwait/alignmentrs/blob/master/LICENSE)

