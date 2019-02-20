use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};

use crate::record::{Record, BaseRecord};

#[pyclass(subclass)]
#[derive(Clone)]
pub struct BaseAlignment {

    pub records: Vec<BaseRecord>,

    pub ncols: i32,
    pub chunk_size: i32,

}

#[pymethods]
impl BaseAlignment {
    #[new]
    /// Creates a new BaseAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, records: Vec<&BaseRecord>, chunk_size: i32) -> PyResult<()> {
        // TODO: Check if all vectors have the same size.
        let ncols = 0;  // TODO: Get ncols from first record else 0
        let mut new_records: Vec<BaseRecord> = Vec::new();
        for b in records.into_iter() {
            let mut record: BaseRecord = b.clone();
            match record.set_chunk_size(chunk_size) {
                Err(x) => return Err(x),
                _ => ()
            };
            new_records.push(record);
        }
        obj.init(|_| {
            BaseAlignment { 
                records: new_records,
                ncols: ncols,
                chunk_size
            }
        })
    }
    
    #[getter]
    pub fn nrows(&self) -> PyResult<i32> {
        Ok(self.records.len() as i32)
    }

    #[getter]
    pub fn ncols(&self) -> PyResult<i32> {
        if self.records.len() == 0 {
            return Ok(0)
        }
        Ok(self.ncols)
    }


    // Sequence getters

    /// get_row(row_index, /)
    /// --
    /// 
    /// Returns a new Record object containing the id, description,
    /// and sequence at the specified index.
    pub fn get_record(&self, row: i32) -> PyResult<BaseRecord> {
        match self.get_records(vec![row]) {
            Ok(vals) => Ok(vals[0].clone()),
            Err(x) => Err(x)
        }
    }

    /// get_rows(row_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified by a list of indices.
    pub fn get_records(&self, rows: Vec<i32>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let mut records: Vec<BaseRecord> = Vec::new();
        for row in rows.into_iter().map(|x| x as usize) {
            check_row_index(self, row)?;
            records.push(self.records[row].clone())
        }
        Ok(records)
    }

    /// get_rows_by_name(names, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified using a list of names/ids.
    pub fn get_records_by_name(&self, names: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = match self.row_names_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_records(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// get_rows_by_prefix(prefixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of prefixes.
    pub fn get_records_by_prefix(&self, prefixes: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = match self.row_prefix_to_indices(prefixes) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_records(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// get_rows_by_suffix(suffixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of suffixes.
    pub fn get_records_by_suffix(&self, suffixes: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = match self.row_suffix_to_indices(suffixes) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_records(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    pub fn get_row(&self, row: i32) -> PyResult<Vec<String>> {
        match self.get_rows(vec![row]) {
            Ok(vals) => Ok(vals[0].clone()),
            Err(x) => Err(x)
        }
    }

    pub fn get_rows(&self, rows: Vec<i32>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let mut sequences: Vec<Vec<String>> = Vec::new();
        for row in rows.into_iter().map(|x| x as usize) {
            check_row_index(self, row)?;
            sequences.push(self.records[row].sequence.clone())
        }
        Ok(sequences)
    }

    pub fn get_rows_by_name(&self, names: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = match self.row_names_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }
    pub fn get_rows_by_prefix(&self, prefixes: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = match self.row_prefix_to_indices(prefixes) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }
    pub fn get_rows_by_suffix(&self, suffixes: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = match self.row_suffix_to_indices(suffixes) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        match self.get_rows(rows) {
            Ok(x) => Ok(x),
            Err(x) => Err(x)
        }
    }

    /// get_col(col_index, /)
    /// --
    /// 
    /// Returns the specified alignment column as a new Record object.
    pub fn get_col(&self, i: i32) -> PyResult<Vec<String>> {
        match self.get_cols(vec![i]) {
            Ok(vals) => Ok(vals[0].clone()),
            Err(x) => Err(x)
        }
    }

    /// get_cols(col_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the specific
    /// alignment columns based on a list of indices. 
    pub fn get_cols(&self, cols: Vec<i32>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let mut records: Vec<Vec<String>> = Vec::new();
        for i in cols.into_iter().map(|x| x as usize) {
            check_col_index(self, i)?;
            let sequence: Vec<String> = self.records.iter()
                .map(|rec| rec.sequence[i].clone()).collect();
            records.push(sequence);
        }
        Ok(records)
    }

    /// subset(row_indices, column_indices, /)
    /// --
    /// 
    /// Returns the subset of rows and columns in the alignment as a new
    /// RawAlignment.
    fn subset(&self, rows: Vec<i32>, cols: Vec<i32>)
    -> PyResult<BaseAlignment> {
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        let mut records: Vec<BaseRecord> = Vec::new();
        let chunk_size = self.chunk_size;
        let ncols = cols.len();
        for row in rows.into_iter().map(|x| x as usize) {
            let sequence: Vec<String> = self.records[row].sequence.iter().enumerate()
                .filter(|(i, _)| row == *i )
                .map(|(_, item)| item.to_string() )
                .collect();
            // Make sure sequence length == ncols
            if sequence.len() != ncols {
                return Err(exceptions::ValueError::py_err(
                    "Unexpected number of columns"))
            }
            records.push(BaseRecord{
                id: self.records[row].id.to_string(),
                description: self.records[row].description.to_string(),
                sequence,
                chunk_size,
            })
        }
        Ok(BaseAlignment{ records, ncols: ncols as i32, chunk_size})
    }


    // Metadata setters

    /// replace_id(row_index, value, /)
    /// --
    /// 
    /// Replaces the name/identifier of an existing entry in the
    /// multiple sequence alignment.
    fn replace_id(&mut self, row: i32, value: &str) -> PyResult<()> {
        self.replace_ids(vec![row], vec![value])
    }

    /// replace_ids(row_indices, values, /)
    /// --
    /// 
    /// Replace the names/identifiers of many entries simultaneously.
    /// Each name/identifier will replace an existing identifier based on the
    /// corresponding list of row indices.
    fn replace_ids(&mut self, rows: Vec<i32>, values: Vec<&str>)
    -> PyResult<()> {
        check_length_match(&rows, &values)?;
        check_empty_alignment(self)?;
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            check_row_index(self, row)?;
            self.records[row].id = values[i].to_string();
        }
        Ok(())
    }

    /// replace_description(index, value, /)
    /// --
    /// 
    /// Replaces the description of an existing entry in the
    /// multiple sequence alignment.
    fn replace_description(&mut self, row: i32, value: &str)
    -> PyResult<()> {
        self.replace_descriptions(vec![row], vec![value])
    }    

    /// replace_descriptions(indices, values, /)
    /// --
    /// 
    /// Replace the descriptions of many entries simultaneously.
    /// Each description will replace an existing description based on the
    /// corresponding list of row indices.
    fn replace_descriptions(&mut self, rows: Vec<i32>, values: Vec<&str>) 
    -> PyResult<()> {
        check_length_match(&rows, &values)?;
        check_empty_alignment(self)?;
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            check_row_index(self, row)?;
            self.records[row].description = values[i].to_string();
        }
        Ok(())
    }

    // Sequence setters

    /// replace_sequence(index, value, /)
    /// --
    ///
    /// Replaces the sequence of an existing entry in the multiple
    /// sequence alignment.
    fn replace_sequence(&mut self, row: i32, value: &str) -> PyResult<()> {
        self.replace_sequences(vec![row], vec![value])
    }

    /// replace_sequences(indices, values, /)
    /// --
    /// 
    /// Replaces many sequences simulateneously.
    /// Each sequence will replace an existing sequence based on the
    /// corresponding list of row indices.
    fn replace_sequences(&mut self, rows: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        check_length_match(&rows, &values)?;
        check_empty_alignment(self)?;
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            check_row_index(self, row)?;
            check_length_match(
                &values[i].chars().collect(), 
                &self.records[i].sequence.join("").chars().collect())?;
            self.records[i].set_sequence(values[i])?;
        }
        Ok(())
    }

    /// replace_col(coordinate, sequence, /)
    /// --
    /// 
    /// Replaces the all the characters in an alignment column.
    fn replace_col(&mut self, col: i32, value: Vec<&str>) -> PyResult<()> {
        self.replace_cols(vec![col], vec![value])
    }

    /// replace_cols(coordinates, sequences, /)
    /// --
    /// 
    /// Replaces the all the characters in each specified column from a list of
    /// alignment columns.
    fn replace_cols(&mut self, cols: Vec<i32>, values: Vec<Vec<&str>>)
    -> PyResult<()> {
        check_length_match(&cols, &values)?;
        check_empty_alignment(self)?;
        for (i, col) in cols.iter().enumerate() {
            check_length_match(&self.records, &values[i])?;
            for row in 0..self.records.len() {
                self.records[row].sequence[*col as usize] = values[i][row].to_string();
            }
        }
        Ok(())
    }

    // Deleters
    // remove_rows and remove_cols are the main methods for deleting
    // contents of RawAlignment

    /// remove_row(index, /)
    /// --
    /// 
    /// Removes one entry corresponding to the specified row index.
    fn remove_row(&mut self, row: i32) -> PyResult<()> {
        self.remove_rows(vec![row])
    }

    /// remove_rows(indices, /)
    /// --
    /// 
    /// Removes many entries simulatenously based on a
    /// list of row indices.
    fn remove_rows(&mut self, mut rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        rows.sort_unstable();
        rows.reverse();
        for row in rows.iter().map(|x| *x as usize) {
            check_row_index(self, row)?;
            self.records.remove(row);
        }
        Ok(())
    }

    /// remove_col(index, /)
    /// --
    /// 
    /// Removes one column corresponding to the specified column index.
    fn remove_col(&mut self, col: i32) -> PyResult<()> {
        self.remove_cols(vec![col])
    }

    /// remove_cols(indices, /)
    /// --
    /// 
    /// Removes many alignment columns simulatenously based on a
    /// list of column indices.
    fn remove_cols(&mut self, mut cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        cols.sort_unstable();
        cols.reverse();
        for row in 0..self.records.len() {
            for col in cols.iter() {
                self.records[row].sequence.remove(*col as usize);
            }
        }
        // Update ncols
        Ok(())
    }

    /// retain_row(index, /)
    /// --
    /// 
    /// Keeps one entry according to the given row index and
    /// removes everything else.
    fn retain_row(&mut self, id: i32) -> PyResult<()> {
        self.retain_rows(vec![id])
    }

    /// retain_rows(indices, /)
    /// 
    /// Keep entries at the specified row indices, and removes
    /// everything else.
    fn retain_rows(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_rows(rows)
    }

    /// retain_col(index, /)
    /// --
    /// 
    /// Keeps one alignment column according to the given column index and
    /// removes everything else.
    fn retain_col(&mut self, col: i32) -> PyResult<()> {
        self.retain_cols(vec![col])
    }

    /// retain_cols(indices, /)
    /// 
    /// Keep  alignment columns at the specified column indices and
    /// removes everything else.
    fn retain_cols(&mut self, cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let cols: Vec<i32> = (0..self.ncols()?)
            .filter(|i| !cols.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();        
        self.remove_cols(cols)
    }

    // The following are extensions of remove_rows and retain_rows
    // that uses sample IDs instead of row indices to reference samples.
    // These are convenience functions that simply do a lookup on the
    // ids vector to get the row ids to use with the remove_row method.

    /// remove_rows_by_name(names, /)
    /// 
    /// Removes samples matching the given sample ID's inplace.
    fn remove_rows_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_names_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        self.remove_rows(rows)
    }

    /// remove_rows_by_prefix(prefixes, /)
    /// 
    /// Removes samples matching at least one of the given prefixes inplace.
    fn remove_rows_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_prefix_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        self.remove_rows(rows)
    }

    /// remove_rows_by_suffix(suffixes, /)
    /// 
    /// Removes samples matching at least one of the given suffixes inplace.
    fn remove_rows_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_suffix_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        self.remove_rows(rows)
    }

    /// retain_rows_by_name(names, /)
    /// 
    /// Keep samples matching the given sample ID's and remove
    /// non-matching samples inplace.
    fn retain_rows_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_names_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_rows(rows)
    }

    /// retain_rows_by_prefix(prefixes, /)
    /// 
    /// Keep samples matching at least one of the given prefixes and remove
    /// non-matching samples inplace.
    fn retain_rows_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_prefix_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_rows(rows)
    }

    /// retain_rows_by_suffix(suffixes, /)
    /// 
    /// Keep samples matching at least one of the given suffixes and remove
    /// non-matching samples inplace.
    fn retain_rows_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = match self.row_suffix_to_indices(names) {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_rows(rows)
    }

    /// insert_row(position, id, description, sequence, /)
    /// --
    /// 
    /// Inserts one entry into the multiple sequence alignment
    /// at the specified position.
    fn insert_row(&mut self, row: i32, value: &BaseRecord) -> PyResult<()> {
        self.insert_rows(vec![row], vec![value])
    }

    /// insert_rows(position, ids, descriptions, sequences, /)
    /// 
    /// Inserts one or more samples at the specified position.
    fn insert_rows(&mut self, rows: Vec<i32>, values: Vec<&BaseRecord>) -> PyResult<()> {
        for (i, value) in values.into_iter().enumerate() {
            let row = rows[i];
            if self.nrows()? < i as i32 {
                return Err(exceptions::IndexError::py_err(
                    "Row index out of range."))
            }
            self.records.insert(i, value.clone());
        }
        Ok(())
    }

    /// append_row(id, description, sequence, /)
    /// --
    /// 
    /// Appends one entry at the end of the multiple sequence alignment.
    fn append_row(&mut self, value: &BaseRecord) -> PyResult<()> {
        self.append_rows(vec![value])
    }

    /// append_rows(ids, descriptions, sequences, /)
    /// 
    /// Appends one or more samples at the end of the list.
    fn append_rows(&mut self, values: Vec<&BaseRecord>) -> PyResult<()> {
        for value in values.into_iter() {
            self.records.push(value.clone());
        }
        Ok(())
    }

    /// reorder_rows(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    fn reorder_rows(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        check_length_match(&rows, &self.records)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
            let mut records: Vec<BaseRecord> = Vec::with_capacity(rows.len());
            for row in rows.iter() {
                records.push(self.records[(*row) as usize].clone())
            }
            self.records = records;
        }
        Ok(())
    }

    /// reorder_cols(ids, /)
    /// --
    /// 
    /// Reorders the alignment columns inplace based on a list of current
    /// column indices.
    fn reorder_cols(&mut self, cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        check_length_match(&cols, &self.records[0].sequence)?;
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
            for i in 0..self.records.len() {
                let sequence = (0..self.records[i].sequence.len())
                    .map(|i| self.records[i].sequence[i].clone())
                    .collect();
                self.records[i].sequence = sequence;
            }
        }
        Ok(())
    }

    // Methods to convert row names to row indices

    /// row_names_to_ids(names, /)
    /// --
    /// 
    /// Converts a list of names to corresponding row indices.
    pub fn row_names_to_indices(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        check_empty_alignment(self)?;
        let mut indices: Vec<i32> = Vec::new();
        for name in names.into_iter() {
            match self.records.iter().position(|b| b.id == name) {
                Some(i) => {
                    indices.push(i as i32);
                },
                None => {
                    return Err(exceptions::ValueError::py_err(
                        format!("sample id {} not found", name)))
                }
            }
        }
        Ok(indices)
    }

    /// row_prefix_to_ids(prefixes, /)
    /// --
    /// 
    /// Matches a list of prefixes to entry names/identifiers and returns
    /// the row indices of matched.
    pub fn row_prefix_to_indices(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        check_empty_alignment(self)?;
        let mut indices: Vec<i32> = Vec::new();
        let mut matches: Vec<&str> = Vec::new();
        for name in names.iter() {
            for (i, b) in self.records.iter().enumerate() {
                if b.id.starts_with(name) && !matches.contains(name) {
                    indices.push(i as i32);
                    matches.push(&b.id);
                }
            }
        }
        Ok(indices)
    }

    /// row_suffix_to_ids(suffixes, /)
    /// --
    /// 
    /// Matches a list of suffixes to entry names/identifiers and returns
    /// the row indices of matched.
    pub fn row_suffix_to_indices(&self, names: Vec<&str>) -> PyResult<Vec<i32>> {
        check_empty_alignment(self)?;
        let mut indices: Vec<i32> = Vec::new();
        let mut matches: Vec<&str> = Vec::new();
        for name in names.iter() {
            for (i, b) in self.records.iter().enumerate() {
                if b.id.ends_with(name) && !matches.contains(name) {
                    indices.push(i as i32);
                    matches.push(&b.id);
                }
            }
        }
        Ok(indices)
    }
    
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseAlignment {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("BaseAlignment(nrows={nrows}, ncols={ncols})",
                   nrows=self.nrows()?, ncols=self.nrows()?))
    }

    fn __str__(&self) -> PyResult<String> {
        if self.nrows()? == 0 {
            return Ok(String::new())
        }
        let mut fasta_strings: Vec<String> = Vec::new();
        for record in self.records.iter() {
            fasta_strings.push(record.__str__()?);
        }
        Ok(fasta_strings.join("\n"))
    }

    // Determines the "truthyness" of the object
    fn __bool__(&self) -> PyResult<bool> {
        if self.ncols()? == 0 {
            return Ok(false)
        }
        Ok(true)
    }
}

// #[pyproto]
// impl PyGCProtocol for BaseAlignment {
//     fn __traverse__(&self, visit: PyVisit) -> Result<(), PyTraverseError> {
//         if self.records.len() > 0 {
//             for obj_ref in self.records.iter() {
//                 visit.call(obj_ref)?
//             }
//         }
//         Ok(())
//     }

//     fn __clear__(&mut self) {
//         if let Some(obj) = self.obj.take() {
//           // Release reference, this decrements ref counter.
//           self.py().release(obj);
//         }
//     }
// }


pub fn check_empty_alignment(aln_ptr: &BaseAlignment) -> PyResult<()> {
    if aln_ptr.nrows()? == 0 {
        return Err(exceptions::ValueError::py_err(
            "Cannot get rows from an empty alignment."))
    }
    Ok(())
}


pub fn check_row_index(aln_ptr: &BaseAlignment, i: usize) -> PyResult<()> {
    if aln_ptr.nrows()? <= i as i32 {
        return Err(exceptions::IndexError::py_err(
            "Row index out of range."))
    }
    Ok(())
}

pub fn check_col_index(aln_ptr: &BaseAlignment, i: usize) -> PyResult<()> {
    if aln_ptr.ncols()? <= i as i32 {
        return Err(exceptions::IndexError::py_err(
            "Column index out of range."))
    }
    Ok(())
}


pub fn check_length_match<T, U>(v1: &Vec<T>, v2: &Vec<U>) -> PyResult<()> {
    if v1.len() != v2.len() {
        return Err(exceptions::ValueError::py_err(
            "Length mismatch."))
    }
    Ok(())
}





























#[pyclass(subclass)]
#[derive(Clone)]
/// RawAlignment(ids, descriptions, sequences)
/// --
/// 
/// RawAlignment represents a multiple sequence alignment.
pub struct RawAlignment {

    #[prop(get)]
    pub ids: Vec<String>,

    #[prop(get)]
    pub descriptions: Vec<String>,
    
    #[prop(get)]
    pub sequences: Vec<String>,  // TODO: Try a Vec<Vec<char>> instead.

}

#[pymethods]
impl RawAlignment {
    #[new]
    /// Creates a new RawAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, ids: Vec<&str>, descriptions: Vec<&str>,
    sequences: Vec<&str>) -> PyResult<()> {
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err(
                "id, description, and sequence lists must have the same length"))
        }
        obj.init(|_| {
            RawAlignment { 
                ids: ids.iter().map(|s| s.to_string()).collect(), 
                descriptions: descriptions.iter().map(|s| s.to_string()).collect(), 
                sequences: sequences.iter().map(|s| s.to_string()).collect() }
        })
    }

    // Sequence getters

    /// get_row(row_index, /)
    /// --
    /// 
    /// Returns a new Record object containing the id, description,
    /// and sequence at the specified index.
    fn get_row(&self, i: usize) -> PyResult<Record> {
        // TODO: Change this to Record to generalize Sample and Marker records
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
        } else if self._nrows() <= i {
            return Err(exceptions::IndexError::py_err(
                "row index out of range"))
        }
        Ok(Record {
            id: self.ids[i].to_string(),
            description: self.descriptions[i].to_string(),
            sequence: self.sequences[i].to_string(),
        })
    }

    /// get_rows(row_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified by a list of indices.
    fn get_rows(&self, ids: Vec<i32>) -> PyResult<Vec<Record>> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
        }
        let mut records: Vec<Record> = Vec::new();
        for i in ids.iter().map(|x| *x as usize) {
            if self._nrows() <= i {
                return Err(exceptions::IndexError::py_err(
                    "row index out of range"))
            } 
            records.push(Record{ 
                id: self.ids[i].to_string(),
                description: self.descriptions[i].to_string(),
                sequence: self.sequences[i].to_string(),
            })
        }
        Ok(records)
    }

    /// get_rows_by_name(names, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified using a list of names/ids.
    fn get_rows_by_name(&self, names: Vec<&str>) -> PyResult<Vec<Record>> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
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

    /// get_rows_by_prefix(prefixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of prefixes.
    fn get_rows_by_prefix(&self, names: Vec<&str>) -> PyResult<Vec<Record>> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
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

    /// get_rows_by_suffix(suffixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of suffixes.
    fn get_rows_by_suffix(&self, names: Vec<&str>) -> PyResult<Vec<Record>> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
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

    /// get_col(col_index, /)
    /// --
    /// 
    /// Returns the specified alignment column as a new Record object.
    fn get_col(&self, i: usize) -> PyResult<Record> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
        }
        if i >= self.sequences[0].chars().count() {
            return Err(exceptions::IndexError::py_err(
                "column index out of range"))
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

    /// get_cols(col_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the specific
    /// alignment columns based on a list of indices. 
    fn get_cols(&self, cols: Vec<i32>) -> PyResult<Vec<Record>> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err(
                "alignment has no sequences"))
        }
        let length = self._ncols();
        let mut records: Vec<Record> = Vec::new();
        for col in cols.into_iter().map(|x| x as usize) {
            if col >= length {
                return Err(exceptions::ValueError::py_err(
                    "column index out of range"))
            }
            let sequence: String = self.sequences.iter().map(|x| {
                let char_vec: Vec<char> = x.chars().collect();
                char_vec[col]
            }).collect();
            records.push(Record{
                id: format!("{}", col),
                description: String::new(),
                sequence,
            });
        }
        Ok(records)
    }

    /// subset(row_indices, column_indices, /)
    /// --
    /// 
    /// Returns the subset of rows and columns in the alignment as a new
    /// RawAlignment.
    fn subset(&self, ids: Vec<i32>, sites: Vec<i32>)
    -> PyResult<RawAlignment> {
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
        Ok(RawAlignment {
            ids: new_ids,
            descriptions: new_descriptions,
            sequences: new_sequences,
        })
    }

    // Metadata setters

    /// replace_id(row_index, value, /)
    /// --
    /// 
    /// Replaces the name/identifier of an existing entry in the
    /// multiple sequence alignment.
    fn replace_id(&mut self, i: i32, value: &str) -> PyResult<()> {
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![value];
        self.replace_ids(ids, values)
    }

    /// replace_ids(row_indices, values, /)
    /// --
    /// 
    /// Replace the names/identifiers of many entries simultaneously.
    /// Each name/identifier will replace an existing identifier based on the
    /// corresponding list of row indices.
    fn replace_ids(&mut self, ids: Vec<i32>, values: Vec<&str>)
    -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err(
                "index and id lists must have the same length"))
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

    /// replace_description(index, value, /)
    /// --
    /// 
    /// Replaces the description of an existing entry in the
    /// multiple sequence alignment.
    fn replace_description(&mut self, i: i32, description: &str)
    -> PyResult<()> {
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![description];
        self.replace_descriptions(ids, values)
    }    

    /// replace_descriptions(indices, values, /)
    /// --
    /// 
    /// Replace the descriptions of many entries simultaneously.
    /// Each description will replace an existing description based on the
    /// corresponding list of row indices.
    fn replace_descriptions(&mut self, ids: Vec<i32>, values: Vec<&str>) 
    -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err(
                "index and description lists must have the same length"))
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

    /// replace_sequence(index, value, /)
    /// --
    ///
    /// Replaces the sequence of an existing entry in the multiple
    /// sequence alignment.
    fn replace_sequence(&mut self, i: i32, seqeunce: &str) -> PyResult<()> {
        let ids: Vec<i32> = vec![i];
        let values: Vec<&str> = vec![seqeunce];
        self.replace_sequences(ids, values)
    }

    /// replace_sequences(indices, values, /)
    /// --
    /// 
    /// Replaces many sequences simulateneously.
    /// Each sequence will replace an existing sequence based on the
    /// corresponding list of row indices.
    fn replace_sequences(&mut self, ids: Vec<i32>, values: Vec<&str>) -> PyResult<()> {
        if ids.len() != values.len() {
            return Err(exceptions::ValueError::py_err(
                "index and sequence lists must have the same length"))
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

    /// replace_col(coordinate, sequence, /)
    /// --
    /// 
    /// Replaces the all the characters in an alignment column.
    fn replace_col(&mut self, col: i32, sequence: &str) -> PyResult<()> {
        if self.ids.len() != sequence.len() {
            return Err(exceptions::ValueError::py_err(
                "sequence and the number of entries in the alignment \
                must have the same length"))
        }
        let col = col as usize;
        let replace_chars_vec: Vec<char> = sequence.chars().collect();
        for row in 0..self.sequences.len() {
            let sequence_vec: String = self.sequences[row].char_indices()
                .map(|(i, y)| {
                    if col == i {
                        replace_chars_vec[row]
                    } else {
                        y
                    }
                }).collect();
            self.sequences[row] = sequence_vec;
        }
        Ok(())
    }

    /// replace_cols(coordinates, sequences, /)
    /// --
    /// 
    /// Replaces the all the characters in each specified column from a list of
    /// alignment columns.
    fn replace_cols(&mut self, cols: Vec<i32>, values: Vec<&str>)
    -> PyResult<()> {
        if cols.len() != values.len() {
            return Err(exceptions::ValueError::py_err(
                "coordinate and sequence lists must have the same length"))
        }
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        // TODO: Get max col value and check if it is less than num of columns
        let replace_chars_vec: Vec<Vec<char>> = values.iter().map(|v| {
            v.chars().collect()
        }).collect();
        let mut cols = cols;
        cols.sort_unstable();
        for row in 0..self.sequences.len() {
            let mut cnt = 0;
            let sequence_vec: String = self.sequences[row].char_indices()
                .map(|(col,y)| {
                    if cols.contains(&(col as i32)) {
                        cnt += 1;
                        replace_chars_vec[cnt-1][row]
                    } else {
                        y
                    }
                }).collect();
            self.sequences[row] = sequence_vec;
        }        
        Ok(())
    }


    // Deleters
    // remove_rows and remove_cols are the main methods for deleting
    // contents of RawAlignment

    /// remove_row(index, /)
    /// --
    /// 
    /// Removes one entry corresponding to the specified row index.
    fn remove_row(&mut self, id: i32) -> PyResult<()> {
        self.remove_rows(vec![id])
    }

    /// remove_rows(indices, /)
    /// --
    /// 
    /// Removes many entries simulatenously based on a
    /// list of row indices.
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

    /// remove_col(index, /)
    /// --
    /// 
    /// Removes one column corresponding to the specified column index.
    fn remove_col(&mut self, id: i32) -> PyResult<()> {
        self.remove_cols(vec![id])
    }

    /// remove_cols(indices, /)
    /// --
    /// 
    /// Removes many alignment columns simulatenously based on a
    /// list of column indices.
    fn remove_cols(&mut self, mut ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        ids.sort_unstable();
        ids.reverse();
        for i in 0..self.sequences.len() {
            let mut sequence_chars: Vec<char> = self.sequences[i].chars()
                .collect();
            for i in ids.iter().map(|x| *x as usize) {
                if i >= sequence_chars.len() {
                    return Err(exceptions::IndexError::py_err("site index out of range"))
                }
                sequence_chars.remove(i);
            }
            self.sequences[i] = sequence_chars.iter().collect();
        }
        Ok(())
    }

    /// retain_row(index, /)
    /// --
    /// 
    /// Keeps one entry according to the given row index and
    /// removes everything else.
    fn retain_row(&mut self, id: i32) -> PyResult<()> {
        self.retain_rows(vec![id])
    }

    /// retain_rows(indices, /)
    /// 
    /// Keep entries at the specified row indices, and removes
    /// everything else.
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

    /// retain_col(index, /)
    /// --
    /// 
    /// Keeps one alignment column according to the given column index and
    /// removes everything else.
    fn retain_col(&mut self, id: i32) -> PyResult<()> {
        self.retain_cols(vec![id])
    }

    /// retain_cols(indices, /)
    /// 
    /// Keep  alignment columns at the specified column indices and
    /// removes everything else.
    fn retain_cols(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        let mut remove_ids: Vec<i32> = Vec::new();
        for i in 0..self.sequences[0].chars().count() {
            if !ids.contains(&(i as i32)) {
                remove_ids.push(i as i32);
            }
        }
        match self.remove_cols(remove_ids) {
            Err(x) => Err(x),
            Ok(x) => Ok(x)
        }
    }

    // The following are extensions of remove_rows and retain_rows
    // that uses sample IDs instead of row indices to reference samples.
    // These are convenience functions that simply do a lookup on the
    // ids vector to get the row ids to use with the remove_row method.

    /// remove_rows_by_name(names, /)
    /// 
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

    /// remove_rows_by_prefix(prefixes, /)
    /// 
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

    /// remove_rows_by_suffix(suffixes, /)
    /// 
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

    /// retain_rows_by_name(names, /)
    /// 
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

    /// retain_rows_by_prefix(prefixes, /)
    /// 
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

    /// retain_rows_by_suffix(suffixes, /)
    /// 
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

    /// insert_row(position, id, description, sequence, /)
    /// --
    /// 
    /// Inserts one entry into the multiple sequence alignment
    /// at the specified position.
    fn insert_row(&mut self, i: i32, id: &str, description: &str, sequence: &str) -> PyResult<()> {
        self.insert_rows(i, vec![id], vec![description], vec![sequence])
    }

    /// insert_rows(position, ids, descriptions, sequences, /)
    /// 
    /// Inserts one or more samples at the specified position.
    fn insert_rows(&mut self, i: i32, ids: Vec<&str>, descriptions: Vec<&str>,
                   sequences: Vec<&str>) -> PyResult<()> {
        let i = i as usize;
        if self._nrows() <= i {
            return Err(exceptions::IndexError::py_err("sample index out of range"))
        }
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err(
                "id, description, and sequence lists must have the same length"))
        }
        for offset in 0..sequences.len() {
            let seq_len = sequences[offset].chars().count();
            if self._nrows() > 0 && self._ncols() != seq_len {
                return Err(exceptions::ValueError::py_err(
                    format!("sequence length does not match the alignment length: {} != {}",
                    seq_len, self._ncols())))
            }
            self.ids.insert(i + offset, ids[offset].to_string());
            self.descriptions.insert(i + offset, descriptions[offset].to_string());
            self.sequences.insert(i + offset, sequences[offset].to_string());
        }
        Ok(())
    }

    /// append_row(id, description, sequence, /)
    /// --
    /// 
    /// Appends one entry at the end of the multiple sequence alignment.
    fn append_row(&mut self, id: &str, description: &str, sequence: &str) -> PyResult<()> {
        self.append_rows(vec![id], vec![description], vec![sequence])
    }

    /// append_rows(ids, descriptions, sequences, /)
    /// 
    /// Appends one or more samples at the end of the list.
    fn append_rows(&mut self, ids: Vec<&str>, descriptions: Vec<&str>,
                   sequences: Vec<&str>) -> PyResult<()> {
        if (ids.len() != descriptions.len()) ||
           (ids.len() != sequences.len()) {
            return Err(exceptions::ValueError::py_err(
                "id, description, and sequence lists must have the same length"))
        }
        for offset in 0..sequences.len() {
            let seq_len = sequences[offset].chars().count();
            if self._nrows() > 0 && self._ncols() != seq_len {
                return Err(exceptions::ValueError::py_err(
                    format!("sequence length does not match the alignment length: {} != {}",
                            seq_len, self._ncols())))
            }
            self.ids.push(ids[offset].to_string());
            self.descriptions.push(descriptions[offset].to_string());
            self.sequences.push(sequences[offset].to_string());
        }
        Ok(())
    }

    /// reorder_rows(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    fn reorder_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if let Some(max) = ids.iter().max() {
            let length = self.ids.len();
            if ids.len() != length {
                return Err(exceptions::ValueError::py_err(format!(
                    "list length does not match the number of samples: {} != {}",
                    ids.len(), length)))
            }
            if *max >= length as i32 {
                return Err(exceptions::IndexError::py_err(
                    format!("index out of range: {}", max)))
            }
            let mut new_ids: Vec<String> = Vec::with_capacity(length);
            let mut new_descriptions: Vec<String> = Vec::with_capacity(length);
            let mut new_sequences: Vec<String> = Vec::with_capacity(length);
            for (i, id) in ids.iter().enumerate() {
                let id = *id as usize;
                new_ids[i].clone_from(&self.ids[id]);
                new_descriptions[i].clone_from(&self.descriptions[id]);
                new_sequences[i].clone_from(&self.sequences[id]);
            }
            self.ids = new_ids;
            self.descriptions = new_descriptions;
            self.sequences = new_sequences;
        }
        Ok(())
    }

    /// reorder_cols(ids, /)
    /// --
    /// 
    /// Reorders the alignment columns inplace based on a list of current
    /// column indices.
    fn reorder_cols(&mut self, ids: Vec<i32>) -> PyResult<()> {
        if self._nrows() == 0 {
            return Err(exceptions::ValueError::py_err("alignment has no sequences"))
        }
        // TODO: Get max of ids and check if it is less than the number of columns
        for i in 0..self.sequences.len() {
            let sequence_vec: Vec<char> = self.sequences[i].chars().collect();
            let new_sequence: String = ids.iter().map(|j| {
                    sequence_vec[*j as usize]
            }).collect();
            self.sequences[i] = new_sequence;
        }
        Ok(())
    }

    /// row_names_to_ids(names, /)
    /// --
    /// 
    /// Converts a list of names to corresponding row indices.
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
                    return Err(exceptions::ValueError::py_err(
                        format!("sample id {} not found", name)))
                }
            }
        }
        Ok(ids)
    }

    /// row_prefix_to_ids(prefixes, /)
    /// --
    /// 
    /// Matches a list of prefixes to entry names/identifiers and returns
    /// the row indices of matched.
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

    /// row_suffix_to_ids(suffixes, /)
    /// --
    /// 
    /// Matches a list of suffixes to entry names/identifiers and returns
    /// the row indices of matched.
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

    /// transpose()
    /// --
    /// 
    /// Transposes the alignment making columns rows and rows columns.
    fn transpose(&self) -> PyResult<RawAlignment> {
        if self._nrows() == 0 {
            let new_ids: Vec<String> = Vec::new();
            let new_descriptions: Vec<String> = Vec::new();
            let new_sequences: Vec<String> = Vec::new();
            return Ok(RawAlignment {
                ids: new_ids,
                descriptions: new_descriptions,
                sequences: new_sequences,
            })
        }
        // Create a matrix representation of the current values
        let old_sequence_matrix: Vec<Vec<char>> = self.sequences.iter()
                                    .map(|x| x.chars().collect()).collect();

        let length = self.sequences[0].chars().count();
        let mut new_ids: Vec<String> = vec![String::new(); length];
        let new_descriptions: Vec<String> = vec![String::new(); length];
        let mut new_sequences: Vec<String> = vec![String::new(); length];
        // Transpose values
        for i in 0..length {
            new_ids[i] = format!("{}", i);
            new_sequences[i] = old_sequence_matrix.iter()
                                .map(|x| x[i]).collect();
        }
        Ok(RawAlignment {
            ids: new_ids,
            descriptions: new_descriptions,
            sequences: new_sequences,
        })
    }

    /// concat(aln_list, /)
    /// --
    /// 
    /// Concatenates a list of alignments to the current alignment side-by-side.
    /// Assumes that the number of rows and the order of rows matches the
    /// current alignment.
    fn concat(&self, aln_list: Vec<&RawAlignment>) -> PyResult<RawAlignment> {
        let ids = self.ids.clone();
        let descriptions = self.descriptions.clone();
        let sequence_len = self._nrows();
        let mut sequences: Vec<String> = vec![String::new(); sequence_len];
        for i in 0..sequence_len {
            sequences[i].push_str(&self.sequences[i]);
            for aln in aln_list.iter() {
                let aln_len = aln.sequences.len();
                if aln_len != sequence_len {
                    return Err(exceptions::ValueError::py_err(
                        format!("cannot concatenate alignments with \
                                 unequal number of samples: {} != {}",
                                 sequence_len, aln_len)))
                }
                sequences[i].push_str(&aln.sequences[i]);
            }
        }
        Ok(RawAlignment {ids, descriptions, sequences})
    }

    /// copy()
    /// --
    /// 
    /// Creates a deep copy of itself.
    fn copy(&self) -> PyResult<RawAlignment> {
        let ids: Vec<String> = self.ids.clone();
        let descriptions: Vec<String> = self.descriptions.clone();
        let sequences: Vec<String> = self.sequences.clone();
        Ok(RawAlignment { ids, descriptions, sequences })
    }

    /// is_row_similar(other, /)
    /// 
    /// Checks if the other alignment has the same number of rows.
    fn is_row_similar(&self, other: &RawAlignment) -> PyResult<bool> {
        if self._nrows() != other._nrows() {
            return Ok(false)
        } else if self.ids != other.ids {
            return Ok(false)
        }
        Ok(true)
    }

    /// is_col_similar(other, /)
    /// 
    /// Checks if the other alignment has the same number of columns.
    fn is_col_similar(&self, other: &RawAlignment) -> PyResult<bool> {
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
    fn shape(&self) -> PyResult<(i32, i32)> {
        Ok((self._nrows() as i32, self._ncols() as i32))
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for RawAlignment {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("RawAlignment(nsamples={nsamples}, ncols={ncols})",
                   nsamples=self._nrows(), ncols=self._ncols()))
    }

    fn __str__(&self) -> PyResult<String> {
        if self._nrows() == 0 {
            return Ok(String::new())
        }
        let mut fasta_strings: Vec<String> = Vec::new();
        for i in 0..self._nrows() {
            if self.descriptions[i].chars().count() > 0 {
                fasta_strings.push(
                    format!(">{} {}\n{}", self.ids[i], self.descriptions[i],
                            self.sequences[i]));
            } else {
                fasta_strings.push(format!(">{}\n{}", self.ids[i], self.sequences[i]));
            }
        }
        Ok(fasta_strings.join("\n"))
    }

    // Determines the "truthyness" of the object
    fn __bool__(&self) -> PyResult<bool> {
        if self._ncols() == 0 {
            return Ok(false)
        }
        Ok(true)
    }
}

impl RawAlignment {
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
/// concat_RawAlignments(aln_list, /)
/// --
/// 
/// Concatenates a list of alignments over the site (column) axis.
fn concat_RawAlignments(aln_list: Vec<&RawAlignment>) -> PyResult<RawAlignment> {
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
                return Err(exceptions::ValueError::py_err(
                    format!("cannot concatenate alignments with unequal \
                            number of samples: {} != {}", 
                            sequence_len, aln_len)))
            }
            sequences[i].push_str(&aln.sequences[i]);
        }
    }
    Ok(RawAlignment {ids, descriptions, sequences})
}

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<RawAlignment>()?;
    m.add_class::<BaseAlignment>()?;
    m.add_function(wrap_function!(concat_RawAlignments))?;

    Ok(())
}
