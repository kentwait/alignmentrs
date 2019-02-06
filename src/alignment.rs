use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

use crate::sample::Sample;

#[pyclass(subclass)]
#[derive(Clone)]
/// BaseAlignment(ids, descriptions, sequences)
/// 
/// BaseAlignment represents a multiple sequence alignment.
pub struct BaseAlignment {

    #[prop(get)]
    pub ids: Vec<String>,

    #[prop(get)]
    pub descriptions: Vec<String>,
    
    #[prop(get)]
    pub sequences: Vec<String>,

}

#[pymethods]
impl BaseAlignment {
    #[new]
    /// Creates a new BaseAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        if (ids.len() != descriptions.len()) ||
           (descriptions.len() != sequences.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err("lists must have the same length"))
        }
        obj.init(|_| {
            BaseAlignment { 
                ids: ids.iter().map(|s| s.to_string()).collect(), 
                descriptions: descriptions.iter().map(|s| s.to_string()).collect(), 
                sequences: sequences.iter().map(|s| s.to_string()).collect() }
        })
    }

    // Sequence getters

    /// Returns the sample id, description, and sequence at the given index as
    /// as a Sample object.
    fn get_sample(&self, i: usize) -> PyResult<Sample> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        } else if i >= self.ids.len() {
            return Err(exceptions::ValueError::py_err("sample index out of range"))
        }
        Ok(Sample {
            id: self.ids[i].to_string(),
            description: self.descriptions[i].to_string(),
            sequence: self.sequences[i].to_string(),
        })
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences by index.
    fn get_samples(&self, ids: Vec<i32>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut new_ids: Vec<String> = Vec::new();
        let mut new_descriptions: Vec<String> = Vec::new();
        let mut new_sequences: Vec<String> = Vec::new();
        for i in ids.iter().map(|x| *x as usize) {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            } 
            new_ids.push(self.ids[i].to_string());
            new_descriptions.push(self.descriptions[i].to_string());
            new_sequences.push(self.sequences[i].to_string());
        }
        Ok(BaseAlignment {
            ids: new_ids,
            descriptions: new_descriptions,
            sequences: new_sequences,
        })
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences by ID.
    fn get_samples_by_name(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences match the given list of prefixes.
    fn get_samples_by_prefix(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences match the given list of suffixes.
    fn get_samples_by_suffix(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns the given site as a Sample object. Uses the given site number
    /// as the sample id of Sample.
    fn get_site(&self, i: usize) -> PyResult<Sample> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        if i >= self.sequences[0].chars().count() {
            return Err(exceptions::ValueError::py_err("site index out of range"))
        }
        let mut site_sequence: Vec<String> = Vec::new();
        for s in self.sequences.iter() {
            let seq: Vec<char> = s.chars().collect();
            site_sequence.push(seq[i].to_string())
        }
        Ok(Sample {
            id: format!("{}", i),
            description: String::new(),
            sequence: site_sequence.join(""),
        })
    }

    /// Returns a new BaseAlignment object containing only the specified
    /// sites.
    fn get_sites(&self, sites: Vec<i32>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut new_sequences: Vec<String> = Vec::new();
        for seq in self.sequences.iter() {
            let mut new_sequence: Vec<String> = Vec::new();
            for i in sites.iter().map(|x| *x as usize) {
                if i >= seq.chars().count() {
                    return Err(exceptions::ValueError::py_err("site index out of range"))
                }
                let seq: Vec<char> = seq.chars().collect();
                new_sequence.push(seq[i].to_string());
            }
            new_sequences.push(new_sequence.join(""))
        }
        Ok(BaseAlignment {
            ids: self.ids.to_vec(),
            descriptions:  self.ids.to_vec(),
            sequences: new_sequences,
        })
    }

    /// Returns the subset of samples and sites of the alignment as a new
    /// BaseAlignment.
    fn subset(&self, ids: Vec<i32>, sites: Vec<i32>) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut new_ids: Vec<String> = Vec::new();
        let mut new_descriptions: Vec<String> = Vec::new();
        let mut new_sequences: Vec<String> = Vec::new();
        for i in ids.iter().map(|x| *x as usize) {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            let mut new_sequence: Vec<String> = Vec::new();
            for i in sites.iter().map(|x| *x as usize) {
                let seq: Vec<char> = self.sequences[i].chars().collect();
                new_sequence.push(seq[i].to_string());
            }
            new_ids.push(self.ids[i].to_string());
            new_descriptions.push(self.descriptions[i].to_string());
            new_sequences.push(new_sequence.join(""))
        }
        Ok(BaseAlignment {
            ids: new_ids,
            descriptions: new_descriptions,
            sequences: new_sequences,
        })
    }

    // Metadata setters

    /// Sets the ID of an existing sample.
    fn set_id(&mut self, i: i32, value: &str) -> PyResult<()> {
        let i = i as usize;
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        } else if i >= self.ids.len() {
            return Err(exceptions::ValueError::py_err("sample index out of range"))
        }
        self.ids[i] = value.to_string();
        Ok(())
    }

    /// Sets the description of an existing sample.
    fn set_description(&mut self, i: i32, description: &str) -> PyResult<()> {
        let i = i as usize;
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        } else if i >= self.ids.len() {
            return Err(exceptions::ValueError::py_err("sample index out of range"))
        }
        self.descriptions[i] = description.to_string();
        Ok(())
    }

    /// Sets the sequence of an existing sample
    fn set_sequence(&mut self, i: i32, seqeunce: &str) -> PyResult<()> {
        let i = i as usize;
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        } else if i >= self.ids.len() {
            return Err(exceptions::ValueError::py_err("sample index out of range"))
        }
        let seq_count = seqeunce.chars().count();
        let aln_len = self.sequences[i].chars().count();
        if seq_count != aln_len {
            return Err(exceptions::ValueError::py_err(format!("sequence length ({}) does not match the alignment length ({})", seq_count, aln_len)))
        }
        self.sequences[i] = seqeunce.to_string();
        Ok(())
    }

    /// Sets many sample IDs simulateneously using a list of indices.
    fn set_ids(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            self.ids[i] = values[c].to_string();
        }
        Ok(())
    }

    /// Sets many sample descriptions simulateneously using a list of indices.
    fn set_descriptions(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            self.descriptions[i] = values[c].to_string();
        }
        Ok(())
    }

    // Sequence setters

    /// Sets many sample sequences simulateneously using a list of indices.
    fn set_sequences(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            if values[c].chars().count() != self.sequences[i].chars().count() {
                return Err(exceptions::ValueError::py_err("sequence length is not the same"))
            }
            self.sequences[i] = values[c].to_string();
        }
        Ok(())
    }

    /// Sets many sample sequences simulateneously using a list of sample IDs.
    fn set_sequences_by_name(&mut self, names: Vec<&str>, values: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            if values[c].chars().count() != self.sequences[i].chars().count() {
                return Err(exceptions::ValueError::py_err("sequence length is not the same"))
            }
            self.sequences[i] = values[c].to_string();
        }
        Ok(())
    }

    // Deleters

    /// Removes samples at the given index positions inplace.
    fn remove_samples(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        ids.sort_unstable();
        ids.reverse();
        for i in ids.iter().map(|x| *x as usize) {
            if i >= self.ids.len() {
                return Err(exceptions::ValueError::py_err("sample index out of range"))
            }
            self.ids.remove(i);
            self.descriptions.remove(i);
            self.sequences.remove(i);
        }
        Ok(())
    }

    /// Removes sites at the specified column positions inplace.
    fn remove_sites(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        ids.sort_unstable();
        ids.reverse();
        for sequence in self.sequences.iter_mut() {
            let mut sequence_chars: Vec<char> = sequence.chars().collect();
            for i in ids.iter().map(|x| *x as usize) {
                if i >= sequence_chars.len() {
                    return Err(exceptions::IndexError::py_err("site index out of range"))
                }
                sequence_chars.remove(i);
            }
            let sequence_str: String = sequence_chars.into_iter().collect();
            *sequence = sequence_str;
        }
        Ok(())
    }

    /// Removes samples matching the given sample ID's inplace.
    fn remove_samples_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Removes samples matching at least one of the given prefixes inplace.
    fn remove_samples_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Removes samples matching at least one of the given suffixes inplace.
    fn remove_samples_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_samples(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Removes sites at the specified column positions inplace.
    fn remove_sites(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self.sequences.len() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        ids.sort_unstable();
        ids.reverse();
        for sequence in self.sequences.iter_mut() {
            let mut sequence_chars: Vec<char> = sequence.chars().collect();
            for i in ids.iter().map(|x| *x as usize) {
                if i >= sequence_chars.len() {
                    return Err(exceptions::ValueError::py_err("site index out of range"))
                }
                sequence_chars.remove(i);
            }
            let sequence_str: String = sequence_chars.into_iter().collect();
            *sequence = sequence_str;
        }
        Ok(())
    }

    /// Keep samples at the given index positions, and remove
    /// non-matching samples inplace.
    fn retain_samples(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_samples(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples matching the given sample ID's and remove
    /// non-matching samples inplace.
    fn retain_samples_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_samples(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples matching at least one of the given prefixes and remove
    /// non-matching samples inplace.
    fn retain_samples_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_samples(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples matching at least one of the given suffixes and remove
    /// non-matching samples inplace.
    fn retain_samples_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.sample_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_samples(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples at the specified column positions and remove
    /// other sites inplace.
    fn retain_sites(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.sequences[0].chars().count() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_sites(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    // TODO: Insert and append sites

    /// Inserts one or more samples at the specified position.
    fn insert_samples(&mut self, i: i32, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        let i = i as usize;
        if self._nsamples() <= i {
            return Err(exceptions::IndexError::py_err("sample index out of range"))
        }
        for offset in 0..sequences.len() {
            self.ids.insert(i + offset, ids[offset].to_string());
            self.descriptions.insert(i + offset, descriptions[offset].to_string());
            self.sequences.insert(i + offset, sequences[offset].to_string());
        }
        Ok(())
    }

    /// Appends one or more samples at the end of the list.
    fn append_samples(&mut self, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        for offset in 0..sequences.len() {
            self.ids.push(ids[offset].to_string());
            self.descriptions.push(descriptions[offset].to_string());
            self.sequences.push(sequences[offset].to_string());
        }
        Ok(())
    }

    /// Converts a list of sample names to its corresponding sample indices.
    fn sample_names_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut ids: Vec<i32> = Vec::new();
        for name in names.iter() {
            match self.ids.iter().position(|x| x == name) {
                Some(i) => {
                    ids.push(i as i32);
                },
                None => {
                    return Err(exceptions::ValueError::py_err(format!("sample id {} not found", name)))
                }
            }
        }
        Ok(ids)
    }

    /// Matches a list of sample prefixes to sample names and returns
    /// the indices of matching samples.
    fn sample_prefix_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut ids: Vec<i32> = Vec::new();
        let mut matches: Vec<&str> = Vec::new();
        for name in names.iter() {
            for (i, id) in self.ids.iter().enumerate() {
                if id.starts_with(name) && !matches.contains(name) {
                    ids.push(i as i32);
                    matches.push(id);
                }
            }
        }
        Ok(ids)
    }

    /// Matches a list of sample suffixes to sample names and returns
    /// the indices of matching samples.
    fn sample_suffix_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nsamples() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut ids: Vec<i32> = Vec::new();
        let mut matches: Vec<&str> = Vec::new();
        for name in names.iter() {
            for (i, id) in self.ids.iter().enumerate() {
                if id.ends_with(name) && !matches.contains(name) {
                    ids.push(i as i32);
                    matches.push(id);
                }
            }
        }
        Ok(ids)
    }

    /// Transposes the alignment making columns rows and rows columns.
    fn transpose(&self) -> PyResult<BaseAlignment> {
        if self._nsamples() == 0 {
            let new_ids: Vec<String> = Vec::new();
            let new_descriptions: Vec<String> = Vec::new();
            let new_sequences: Vec<String> = Vec::new();
            return Ok(BaseAlignment {
                ids: new_ids,
                descriptions: new_descriptions,
                sequences: new_sequences,
            })
        }
        let length = self.sequences[0].chars().count();
        let mut new_ids: Vec<String> = Vec::new();
        let new_descriptions: Vec<String> = vec![String::new(); length];
        let mut new_sequences: Vec<String> = Vec::new();
        
        for i in 0..length {
            let mut column = String::new();
            for seq in self.sequences.iter() {
                let sequence_vec: Vec<char> = seq.chars().collect();
                column.push(sequence_vec[i]);
            }
            new_ids.push(format!("{}", i));
            new_sequences.push(column);
        }
        Ok(BaseAlignment {
            ids: new_ids,
            descriptions: new_descriptions,
            sequences: new_sequences,
        })
    }

    /// Concatenates a list of alignments to the current alignment side-by-site
    /// on the site (column) axis.
    fn concat(&self, aln_list: Vec<&BaseAlignment>) -> PyResult<BaseAlignment> {
        let ids = self.ids.clone();
        let descriptions = self.descriptions.clone();
        let sequence_len = self._nsamples();
        let mut sequences: Vec<String> = vec![String::new(); sequence_len];
        for i in 0..sequence_len {
            sequences[i].push_str(&self.sequences[i]);
            for aln in aln_list.iter() {
                let aln_len = aln.sequences.len();
                if aln_len != sequence_len {
                    return Err(exceptions::ValueError::py_err(format!("cannot concatenate alignments with unequal number of samples: {} != {}", sequence_len, aln_len)))
                }
                sequences[i].push_str(&aln.sequences[i]);
            }
        }
        Ok(BaseAlignment {ids, descriptions, sequences})
    }

    // Properties
    #[getter]
    fn nsamples(&self) -> PyResult<i32> {
        Ok(self._nsamples() as i32)
    }

    #[getter]
    fn nsites(&self) -> PyResult<i32> {
        Ok(self._nsites() as i32)
    }

    #[getter]
    fn shape(&self) -> PyResult<(i32, i32)> {
        Ok((self._nsamples() as i32, self._nsites() as i32))
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseAlignment {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("BaseAlignment(nsamples={nsamples}, nsites={nsites})", nsamples=self._nsamples(), nsites=self._nsites()))
    }

    fn __str__(&self) -> PyResult<String> {
        let mut fasta_strings: Vec<String> = Vec::new();
        for i in 0..self._nsamples() {
            if self.descriptions[i].chars().count() > 0 {
                fasta_strings.push(format!(">{} {}\n{}", self.ids[i], self.descriptions[i], self.sequences[i]));
            } else {
                fasta_strings.push(format!(">{}\n{}", self.ids[i], self.sequences[i]));
            }
        }
        Ok(fasta_strings.join("\n"))
    }
}

impl BaseAlignment {
    fn _nsamples(&self) -> usize {
        self.ids.len()
    }

    fn _nsites(&self) -> usize {
        match self.ids.len() {
            x if x == 0 => 0,
            _ => self.sequences[0].chars().count(),
        }
    }
}

#[pyfunction]
/// Concatenates a list of alignments over the site (column) axis.
fn concat_basealignments(aln_list: Vec<&BaseAlignment>) -> PyResult<BaseAlignment> {
    if aln_list.len() == 0 {
        return Err(exceptions::ValueError::py_err("empty list"))
    }
    let ids = aln_list[0].ids.clone();
    let descriptions = aln_list[0].descriptions.clone();
    let sequence_len = aln_list[0].sequences.len();
    let mut sequences: Vec<String> = vec![String::new(); sequence_len];
    for i in 0..sequence_len {
        for aln in aln_list.iter() {
            let aln_len = aln.sequences.len();
            if aln_len != sequence_len {
                return Err(exceptions::ValueError::py_err(format!("cannot concatenate alignments with unequal number of samples: {} != {}", sequence_len, aln_len)))
            }
            sequences[i].push_str(&aln.sequences[i]);
        }
    }
    Ok(BaseAlignment {ids, descriptions, sequences})
}

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(concat_basealignments))?;

    // Add Block class
    m.add_class::<BaseAlignment>()?;

    Ok(())
}