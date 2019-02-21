use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};

use crate::record::BaseRecord;

#[pyclass(subclass)]
#[derive(Clone)]
pub struct BaseAlignment {

    #[prop(get)]
    pub records: Vec<BaseRecord>,

    #[prop(get)]
    pub chunk_size: i32,

}

#[pymethods]
impl BaseAlignment {
    #[new]
    /// Creates a new BaseAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, records: Vec<&BaseRecord>, chunk_size: i32) -> PyResult<()> {
        // TODO: Check if all vectors have the same size.
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
                chunk_size,
            }
        })
    }
    
    #[getter]
    /// int: Returns the number of rows in the BaseAlignment.
    pub fn nrows(&self) -> PyResult<i32> {
        Ok(self.records.len() as i32)
    }

    #[getter]
    /// int: Returns the number of columns in the alignment.
    pub fn ncols(&self) -> PyResult<i32> {
        if self.records.len() == 0 {
            return Ok(0)
        }
        self.records[0].len()
    }

    #[getter]
    /// int: Returns the number of aligned characters (ncols * chunk_size).
    pub fn nchars(&self) -> PyResult<i32> {
        if self.records.len() == 0 {
            return Ok(0)
        }
        self.records[0].str_len()
    }

    #[getter]
    /// list of str: Returns the list of identifiers.
    pub fn ids(&self) -> PyResult<Vec<String>> {
        let values: Vec<String> = self.records.iter()
            .map(|r| r.id.to_string() )
            .collect();
        Ok(values)
    }

    #[getter]
    /// list of str: Returns the list of descriptions.
    pub fn descriptions(&self) -> PyResult<Vec<String>> {
        let values: Vec<String> = self.records.iter()
            .map(|r| r.description.to_string() )
            .collect();
        Ok(values)
    }

    #[getter]
    /// list of str: Returns the list of sequences.
    pub fn sequences(&self) -> PyResult<Vec<String>> {
        let values: Vec<String> = self.records.iter()
            .map(|r| r.sequence.join("") )
            .collect();
        Ok(values)
    }

    #[getter]
    /// list of str: Returns the list of sequences.
    pub fn chunked_sequences(&self) -> PyResult<Vec<Vec<String>>> {
        let values: Vec<Vec<String>> = self.records.iter()
            .map(|r| r.sequence.clone() )
            .collect();
        Ok(values)
    }

    #[setter(chunk_size)]
    pub fn set_chunk_size(&mut self, chunk_size: i32) -> PyResult<()> {
        for i in 0..self.records.len() {
            self.records[i].set_chunk_size(chunk_size)?;
        }
        Ok(())
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

    pub fn get_chunk(&self, i: usize, chunk_size: usize) -> PyResult<Vec<String>> {
        let sequences: Vec<String> = self.records.iter()
                .map(|rec| rec.sequence[i..i+chunk_size].join("")).collect();
        Ok(sequences)
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
        Ok(BaseAlignment{ records, chunk_size })
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

    fn replace_row(&mut self, row: i32, value: &BaseRecord) -> PyResult<()> {
        self.replace_rows(vec![row], vec![value])
    }

    fn replace_rows(&mut self, rows: Vec<i32>, values: Vec<&BaseRecord>) -> PyResult<()> {
        check_length_match(&rows, &values)?;
        check_empty_alignment(self)?;
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            check_row_index(self, row)?;
            // TODO: Make function that checks if records have the same length and string length
            check_length_match_i32(values[i].str_len()?, self.records[i].str_len()?)?;
            check_length_match_i32(values[i].len()?, self.records[i].len()?)?;
            self.records[row] = values[i].clone();
        }
        Ok(())
    }

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
            check_length_match_i32(values[i].len() as i32, self.records[i].str_len()?)?;
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
        self.insert_rows(row, vec![value])
    }

    /// insert_rows(position, ids, descriptions, sequences, /)
    /// 
    /// Inserts one or more samples at the specified position.
    fn insert_rows(&mut self, mut row: i32, values: Vec<&BaseRecord>) -> PyResult<()> {
        if self.nrows()? < row {
            return Err(exceptions::IndexError::py_err(
                "Row index out of range."))
        }
        for value in values.into_iter() {
            self.records.insert(row as usize, value.clone());
            row += 1;
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
                let sequence = cols.iter()
                    .map(|j| self.records[i].sequence[(*j) as usize].to_string())
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
    
    pub fn invert_rows(&self, rows: Vec<i32>) -> PyResult<Vec<i32>> {
        let rows: Vec<i32> = (0..self.nrows()?)
                .filter(|i| !rows.contains(&(*i as i32)) )
                .map(|i| i as i32 )
                .collect();
        Ok(rows)
    }

    pub fn invert_cols(&self, cols: Vec<i32>) -> PyResult<Vec<i32>> {
        let cols: Vec<i32> = (0..self.ncols()?)
                .filter(|i| !cols.contains(&(*i as i32)) )
                .map(|i| i as i32 )
                .collect();
        Ok(cols)
    }

    pub fn copy(&self) -> PyResult<BaseAlignment> {
        Ok(BaseAlignment{ records: self.records.clone(), chunk_size: self.chunk_size })
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseAlignment {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("BaseAlignment(nrows={nrows}, ncols={ncols})",
                   nrows=self.nrows()?, ncols=self.ncols()?))
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

pub fn check_length_match_i32(len1: i32, len2:i32) -> PyResult<()> {
    if len1 != len2 {
        return Err(exceptions::ValueError::py_err(
            "Length mismatch."))
    }
    Ok(())
}

// TODO: Make concat_basealignments

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BaseAlignment>()?;

    Ok(())
}
