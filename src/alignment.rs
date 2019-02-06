use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

use crate::record::Record;

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
    pub sequences: Vec<String>,  // TODO: Try a Vec<Vec<char>> instead.

}

#[pymethods]
impl BaseAlignment {
    #[new]
    /// Creates a new BaseAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err("id, description, and sequence lists must have the same length"))
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
    fn get_row(&self, i: usize) -> PyResult<Record> {  // TODO: Change this to Record to generalize Sample and Marker records
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        } else if self._nrows() <= i {
            return Err(exceptions::IndexError::py_err("sample index out of range"))
        }
        Ok(Record {
            id: self.ids[i].to_string(),
            description: self.descriptions[i].to_string(),
            sequence: self.sequences[i].to_string(),
        })
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences by index.
    fn get_rows(&self, ids: Vec<i32>) -> PyResult<BaseAlignment> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut new_ids: Vec<String> = Vec::new();
        let mut new_descriptions: Vec<String> = Vec::new();
        let mut new_sequences: Vec<String> = Vec::new();
        for i in ids.iter().map(|x| *x as usize) {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
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
    fn get_rows_by_name(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences match the given list of prefixes.
    fn get_rows_by_prefix(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns a new BaseAlignment object containing the specified
    /// sample sequences match the given list of suffixes.
    fn get_rows_by_suffix(&self, names: Vec<&str>) -> PyResult<BaseAlignment> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Returns the given site as a Sample object. Uses the given site number
    /// as the sample id of Sample.
    fn get_site(&self, i: usize) -> PyResult<Record> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        if i >= self.sequences[0].chars().count() {
            return Err(exceptions::IndexError::py_err("site index out of range"))
        }
        let mut site_sequence: Vec<String> = Vec::new();
        for s in self.sequences.iter() {
            let seq: Vec<char> = s.chars().collect();
            site_sequence.push(seq[i].to_string())
        }
        Ok(Record {
            id: format!("{}", i),
            description: String::new(),
            sequence: site_sequence.join(""),
        })
    }

    /// Returns a new BaseAlignment object containing only the specified
    /// sites.
    fn get_sites(&self, sites: Vec<i32>) -> PyResult<BaseAlignment> {
        if self._nrows() == 0 {
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
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        match ids.iter().max() {
            Some(x) if *x as usize >= self._nrows() => {
                return Err(exceptions::IndexError::py_err(
                    "row position cannot be larger than the number of rows"))
            },
            _ => ()
        }
        match sites.iter().max() {
            Some(x) if *x as usize >= self._ncols() => {
                return Err(exceptions::IndexError::py_err(
                    "site position cannot be larger than the number of sites"))
            },
            _ => ()
        }
        let mut new_ids: Vec<String> = Vec::new();
        let mut new_descriptions: Vec<String> = Vec::new();
        let mut new_sequences: Vec<String> = Vec::new();
        for i in ids.iter().map(|x| *x as usize) {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
            }
            let new_sequence: String = self.sequences[i].chars().enumerate()
                                        .filter(|(i, _)| sites.contains(&(*i as i32)))
                                        .map(|(_, x)| x )
                                        .collect();
            new_ids.push(self.ids[i].to_string());
            new_descriptions.push(self.descriptions[i].to_string());
            new_sequences.push(new_sequence)
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
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![value];
        match self.set_ids(ids, values) {
            Err(x) => Err(x),
            _ => Ok(()),
        }
    }

    /// Sets many sample IDs simulateneously using a list of indices.
    fn set_ids(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err("index and id lists must have the same length"))
        }
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
            }
            self.ids[i] = values[c].to_string();
        }
        Ok(())
    }

    /// Sets the description of an existing sample.
    fn set_description(&mut self, i: i32, description: &str) -> PyResult<()> {
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![description];
        match self.set_descriptions(ids, values) {
            Err(x) => Err(x),
            _ => Ok(()),
        }
    }    

    /// Sets many sample descriptions simulateneously using a list of indices.
    fn set_descriptions(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err("index and description lists must have the same length"))
        }
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
            }
            self.descriptions[i] = values[c].to_string();
        }
        Ok(())
    }

    // Sequence setters

    /// Sets the sequence of an existing sample
    fn set_sequence(&mut self, i: i32, seqeunce: &str) -> PyResult<()> {
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![seqeunce];
        match self.set_sequences(ids, values) {
            Err(x) => Err(x),
            _ => Ok(()),
        }
    }

    /// Sets many sample sequences simulateneously using a list of indices.
    fn set_sequences(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err("index and sequence lists must have the same length"))
        }
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
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
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err(format!("number of matched rows is not equal to the length of the given sequence list: {} != {}", ids.len(), values.len())))
        }
        for (c, i) in ids.into_iter().map(|x| x as usize).enumerate() {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
            }
            if values[c].chars().count() != self.sequences[i].chars().count() {
                return Err(exceptions::ValueError::py_err("sequence length is not the same"))
            }
            self.sequences[i] = values[c].to_string();
        }
        Ok(())
    }

    // Deleters
    // remove_rows and remove_sites are the main methods for deleting
    // contents of BaseAlignment

    /// Removes samples at the given index positions inplace.
    /// Index positions are specified by a list of integer ids.
    fn remove_rows(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        ids.sort_unstable();
        ids.reverse();
        for i in ids.iter().map(|x| *x as usize) {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err("sample index out of range"))
            }
            self.ids.remove(i);
            self.descriptions.remove(i);
            self.sequences.remove(i);
        }
        Ok(())
    }

    /// Removes sites at the specified column positions inplace.
    fn remove_sites(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
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

    /// Keep samples at the given index positions, and remove
    /// non-matching samples inplace.
    /// This is the opposite of `remove_rows(ids)`.
    fn retain_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_rows(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples at the specified column positions and remove
    /// other sites inplace.
    /// This is the opposite of `remove_sites(ids)`.
    fn retain_sites(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
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

    // The following are extensions of remove_rows and retain_rows
    // that uses sample IDs instead of row indices to reference samples.
    // These are convenience functions that simply do a lookup on the
    // ids vector to get the row ids to use with the remove_row method.

    /// Removes samples matching the given sample ID's inplace.
    fn remove_rows_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Removes samples matching at least one of the given prefixes inplace.
    fn remove_rows_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Removes samples matching at least one of the given suffixes inplace.
    fn remove_rows_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.remove_rows(ids) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// Keep samples matching the given sample ID's and remove
    /// non-matching samples inplace.
    fn retain_rows_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_names_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_rows(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples matching at least one of the given prefixes and remove
    /// non-matching samples inplace.
    fn retain_rows_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_prefix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_rows(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    /// Keep samples matching at least one of the given suffixes and remove
    /// non-matching samples inplace.
    fn retain_rows_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let ids = match self.row_suffix_to_ids(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.ids.len() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_rows(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    // TODO: Insert and append sites
    // TODO: Insert/append ONE sample using Record object

    /// Inserts one or more samples at the specified position.
    fn insert_rows(&mut self, i: i32, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        let i = i as usize;
        if self._nrows() <= i {
            return Err(exceptions::IndexError::py_err("sample index out of range"))
        }
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err("id, description, and sequence lists must have the same length"))
        }
        for offset in 0..sequences.len() {
            let seq_len = sequences[offset].chars().count();
            if self._nrows() > 0 && self._ncols() != seq_len {
                return Err(exceptions::ValueError::py_err(format!("sequence length does not match the alignment length: {} != {}", seq_len, self._ncols())))
            }
            self.ids.insert(i + offset, ids[offset].to_string());
            self.descriptions.insert(i + offset, descriptions[offset].to_string());
            self.sequences.insert(i + offset, sequences[offset].chars().count().to_string());
        }
        Ok(())
    }

    /// Appends one or more samples at the end of the list.
    fn append_rows(&mut self, ids: Vec<&str>, descriptions: Vec<&str>, sequences: Vec<&str>) -> PyResult<()> {
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err("id, description, and sequence lists must have the same length"))
        }
        for offset in 0..sequences.len() {
            let seq_len = sequences[offset].chars().count();
            if self._nrows() > 0 && self._ncols() != seq_len {
                return Err(exceptions::ValueError::py_err(format!("sequence length does not match the alignment length: {} != {}", seq_len, self._ncols())))
            }
            self.ids.push(ids[offset].to_string());
            self.descriptions.push(descriptions[offset].to_string());
            self.sequences.push(sequences[offset].to_string());
        }
        Ok(())
    }

    /// Converts a list of sample names to its corresponding sample indices.
    fn row_names_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nrows() == 0 {
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
    fn row_prefix_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nrows() == 0 {
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
    fn row_suffix_to_ids(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        if self._nrows() == 0 {
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
        if self._nrows() == 0 {
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
        let mut new_sequences_matrix: Vec<Vec<char>> = vec![
            vec!['#'; self.sequences.len()]; length];
        // Transpose values
        for (i, seq) in self.sequences.iter().enumerate() {
            let sequence_vec: Vec<char> = seq.chars().collect();
            for j in 0..length {
                new_sequences_matrix[j][i] = sequence_vec[i];
            }
        }
        // Finalize
        let mut new_sequences: Vec<String> = Vec::new();
        for j in 0..length {
            new_ids.push(format!("{}", j));
            let seq: String = new_sequences_matrix[j].iter().collect();;
            new_sequences.push(seq);
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
        let sequence_len = self._nrows();
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

    fn is_row_similar(&self, other: &BaseAlignment) -> PyResult<bool> {
        if self._nrows() != other._nrows() {
            return Ok(false)
        } else if self.ids != other.ids {
            return Ok(false)
        }
        Ok(true)
    }

    fn is_col_similar(&self, other: &BaseAlignment) -> PyResult<bool> {
        if self._ncols() != other._ncols() {
            return Ok(false)
        }
        Ok(true)
    }

    // Properties
    #[getter]
    fn nrows(&self) -> PyResult<i32> {
        Ok(self._nrows() as i32)
    }

    #[getter]
    fn ncols(&self) -> PyResult<i32> {
        Ok(self._ncols() as i32)
    }

    #[getter]
    fn nsites(&self) -> PyResult<i32> {
        Ok(self._ncols() as i32)
    }

    #[getter]
    fn shape(&self) -> PyResult<(i32, i32)> {
        Ok((self._nrows() as i32, self._ncols() as i32))
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseAlignment {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("BaseAlignment(nsamples={nsamples}, nsites={nsites})", nsamples=self._nrows(), nsites=self._ncols()))
    }

    fn __str__(&self) -> PyResult<String> {
        if self._nrows() == 0 {
            return Ok(String::new())
        }
        let mut fasta_strings: Vec<String> = Vec::new();
        for i in 0..self._nrows() {
            if self.descriptions[i].chars().count() > 0 {
                fasta_strings.push(format!(">{} {}\n{}", self.ids[i], self.descriptions[i], self.sequences[i]));
            } else {
                fasta_strings.push(format!(">{}\n{}", self.ids[i], self.sequences[i]));
            }
        }
        Ok(fasta_strings.join("\n"))
    }

    // Determines the "truthyness" of the object
    fn __bool__(&self) -> PyResult<bool> {
        if self._nrows() == 0 {
            return Ok(false)
        }
        Ok(true)
    }
}

impl BaseAlignment {
    fn _nrows(&self) -> usize {
        self.ids.len()
    }

    fn _ncols(&self) -> usize {
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