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
        self.records[0].len_str()
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
        self.chunk_size = chunk_size;
        Ok(())
    }


    // Row methods

    /// get_record(row_index, /)
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

    /// get_records(row_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified by a list of indices.
    pub fn get_records(&self, rows: Vec<i32>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut records: Vec<BaseRecord> = Vec::new();
        for row in rows.into_iter().map(|x| x as usize) {
            records.push(self.records[row].clone());
        }
        Ok(records)
    }

    /// get_records_by_name(names, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// specified using a list of names/ids.
    pub fn get_records_by_name(&self, names: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = self.row_names_to_indices(names)?;
        self.get_records(rows)
    }

    /// get_records_by_prefix(prefixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of prefixes.
    pub fn get_records_by_prefix(&self, prefixes: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = self.row_prefix_to_indices(prefixes)?;
        self.get_records(rows)
    }

    /// get_records_by_suffix(suffixes, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the sequences
    /// that match the given list of suffixes.
    pub fn get_records_by_suffix(&self, suffixes: Vec<&str>) -> PyResult<Vec<BaseRecord>> {
        check_empty_alignment(self)?;
        let rows = self.row_suffix_to_indices(suffixes)?;
        self.get_records(rows)
    }

    pub fn get_row(&self, row: i32) -> PyResult<Vec<String>> {
        match self.get_rows(vec![row]) {
            Ok(vals) => Ok(vals[0].clone()),
            Err(x) => Err(x)
        }
    }

    pub fn get_rows(&self, rows: Vec<i32>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut sequences: Vec<Vec<String>> = Vec::new();
        for row in rows.into_iter().map(|x| x as usize) {
            sequences.push(self.records[row].sequence.clone());
        }
        Ok(sequences)
    }

    pub fn get_rows_by_name(&self, names: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = self.row_names_to_indices(names)?;
        self.get_rows(rows)
    }

    pub fn get_rows_by_prefix(&self, prefixes: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = self.row_prefix_to_indices(prefixes)?;
        self.get_rows(rows)
    }

    pub fn get_rows_by_suffix(&self, suffixes: Vec<&str>) -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        let rows = self.row_suffix_to_indices(suffixes)?;
        self.get_rows(rows)
    }

    /// insert_record(position, id, description, sequence, /)
    /// --
    /// 
    /// Inserts one entry into the multiple sequence alignment
    /// at the specified position.
    fn insert_record(&mut self, row: i32, value: &BaseRecord) -> PyResult<()> {
        self.insert_records(row, vec![value])
    }

    /// insert_records(position, ids, descriptions, sequences, /)
    /// 
    /// Inserts one or more samples at the specified position.
    fn insert_records(&mut self, mut row: i32, values: Vec<&BaseRecord>) -> PyResult<()> {
        if self.nrows()? < row {
            return Err(exceptions::IndexError::py_err(
                "row index out of range."))
        }
        for value in values.into_iter() {
            self.records.insert(row as usize, value.clone());
            row += 1;
        }
        Ok(())
    }

    /// append_record(id, description, sequence, /)
    /// --
    /// 
    /// Appends one entry at the end of the multiple sequence alignment.
    fn append_record(&mut self, value: &BaseRecord) -> PyResult<()> {
        self.append_records(vec![value])
    }

    /// append_records(ids, descriptions, sequences, /)
    /// 
    /// Appends one or more samples at the end of the list.
    fn append_records(&mut self, values: Vec<&BaseRecord>) -> PyResult<()> {
        for value in values.into_iter() {
            self.records.push(value.clone());
        }
        Ok(())
    }

    /// remove_record(index, /)
    /// --
    /// 
    /// Removes one entry corresponding to the specified row index.
    fn remove_record(&mut self, row: i32) -> PyResult<()> {
        self.remove_records(vec![row])
    }

    /// remove_records(indices, /)
    /// --
    /// 
    /// Removes many entries simulatenously based on a
    /// list of row indices.
    fn remove_records(&mut self, mut rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        rows.sort_unstable();
        rows.dedup();
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        rows.reverse();
        for row in rows.iter().map(|x| *x as usize) {
            self.records.remove(row);
        }
        Ok(())
    }

    /// remove_records_by_name(names, /)
    /// 
    /// Removes samples matching the given sample ID's inplace.
    fn remove_records_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_names_to_indices(names)?;
        self.remove_records(rows)
    }

    /// remove_records_by_prefix(prefixes, /)
    /// 
    /// Removes samples matching at least one of the given prefixes inplace.
    fn remove_records_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_prefix_to_indices(names)?;
        self.remove_records(rows)
    }

    /// remove_records_by_suffix(suffixes, /)
    /// 
    /// Removes samples matching at least one of the given suffixes inplace.
    fn remove_records_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_suffix_to_indices(names)?;
        self.remove_records(rows)
    }

    /// retain_record(index, /)
    /// --
    /// 
    /// Keeps one entry according to the given row index and
    /// removes everything else.
    fn retain_record(&mut self, id: i32) -> PyResult<()> {
        self.retain_records(vec![id])
    }

    /// retain_records(indices, /)
    /// 
    /// Keep entries at the specified row indices, and removes
    /// everything else.
    fn retain_records(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows: Vec<i32> = self.invert_rows(rows)?;
        self.remove_records(rows)
    }

    fn drain_records(&mut self, mut rows: Vec<i32>) -> PyResult<BaseAlignment> {
        check_empty_alignment(self)?;
        rows.sort_unstable();
        rows.dedup();
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut records: Vec<BaseRecord> = Vec::new();
        rows.reverse();
        for row in rows.into_iter().map(|x| x as usize) {
            records.insert(0, self.records[row].clone());
            self.records.remove(row);
        }
        Ok(BaseAlignment{ records, chunk_size: self.chunk_size}) 
    }

    fn replace_record(&mut self, row: i32, value: &BaseRecord) -> PyResult<()> {
        self.replace_records(vec![row], vec![value])
    }

    fn replace_records(&mut self, rows: Vec<i32>, values: Vec<&BaseRecord>) -> PyResult<()> {
        check_length_match(&rows, &values)?;
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            // TODO: Make function that checks if records have the same length and string length
            check_length_match_i32(values[i].len_str()?, self.records[i].len_str()?)?;
            check_length_match_i32(values[i].len()?, self.records[i].len()?)?;
            self.records[row] = values[i].clone();
        }
        Ok(())
    }

    /// reorder_records(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    fn reorder_records(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_length_match(&rows, &self.records)?;
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut records: Vec<BaseRecord> = Vec::with_capacity(rows.len());
            for row in rows.iter() {
                records.push(self.records[(*row) as usize].clone());
            }
            self.records = records;
        Ok(())
    }



    // Column methods
    

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
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        let mut records: Vec<Vec<String>> = Vec::new();
        for i in cols.into_iter().map(|x| x as usize) {
            let sequence: Vec<String> = self.records.iter()
                .map(|rec| rec.sequence[i].clone())
                .collect();
            records.push(sequence);
        }
        Ok(records)
    }

    pub fn get_chunk(&self, i: usize, chunk_size: usize) -> PyResult<Vec<String>> {
        let sequences: Vec<String> = self.records.iter()
            .map(|rec| {
                let chars: Vec<char> = rec.sequence.join("").chars().collect();
                chars[i..i+chunk_size].into_iter().collect::<String>()
            })
            .collect();
        Ok(sequences)
    }

    // insert
    fn insert_col(&mut self, col: usize, value: Vec<&str>) -> PyResult<()> {
        self.insert_cols(col, vec![value])
    }

    fn insert_cols(&mut self, mut col: usize, values: Vec<Vec<&str>>)
    -> PyResult<()> {
        if values.len() == 0 {
            return Ok(())
        }
        // TODO: Add immediate return to other methods when value length is 0
        check_col_index(self, col)?;
        check_length_match(&values[0], &self.records)?;
        check_chunk_size(self, &values[0])?;
        for value in values.iter() {
            for i in 0..self.records.len() {
                self.records[i].sequence.insert(
                    col as usize, value[i].to_string());
            }
            col += 1;
        }
        Ok(())
    }

    // for value in values.into_iter() {
    //     self.records.insert(row as usize, value.clone());
    //     row += 1;
    // }
    // Ok(())

    // append
    fn append_col(&mut self, value: Vec<&str>) -> PyResult<()> {
        self.append_cols(vec![value])
    }

    fn append_cols(&mut self, values: Vec<Vec<&str>>) -> PyResult<()> {
        if values.len() == 0 {
            return Ok(())
        }
        check_length_match(&values[0], &self.records)?;
        check_chunk_size(self, &values[0])?;
        for i in 0..self.records.len() {
            for value in values.iter() {
                self.records[i].sequence.push(value[i].to_string());
            }
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
        cols.dedup();
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        cols.reverse();
        for row in 0..self.records.len() {
            for col in cols.iter() {
                self.records[row].sequence.remove(*col as usize);
            }
        }
        Ok(())
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

    fn drain_cols(&mut self, cols: Vec<i32>) -> PyResult<BaseAlignment> {
        let mut aln = self.copy()?;
        aln.retain_cols(cols.clone())?;
        self.remove_cols(cols.clone())?;
        Ok(aln)
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
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        for (i, col) in cols.iter().enumerate() {
            check_length_match(&self.records, &values[i])?;
            for row in 0..self.records.len() {
                self.records[row].sequence[*col as usize] = values[i][row].to_string();
            }
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
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        check_length_match(&cols, &self.records[0].sequence)?;
        for i in 0..self.records.len() {
            let sequence = cols.iter()
                .map(|j| self.records[i].sequence[(*j) as usize].to_string())
                .collect();
            self.records[i].sequence = sequence;
        }
        Ok(())
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
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
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
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            self.records[row].description = values[i].to_string();
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
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
            check_length_match_i32(values[i].len() as i32, self.records[i].len_str()?)?;
            self.records[row].set_sequence(values[i])?;
        }
        Ok(())
    }

    /// retain_records_by_name(names, /)
    /// 
    /// Keep samples matching the given sample ID's and remove
    /// non-matching samples inplace.
    fn retain_records_by_name(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_names_to_indices(names)?;
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_records(rows)
    }

    /// retain_records_by_prefix(prefixes, /)
    /// 
    /// Keep samples matching at least one of the given prefixes and remove
    /// non-matching samples inplace.
    fn retain_records_by_prefix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_prefix_to_indices(names)?;
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_records(rows)
    }

    /// retain_records_by_suffix(suffixes, /)
    /// 
    /// Keep samples matching at least one of the given suffixes and remove
    /// non-matching samples inplace.
    fn retain_records_by_suffix(&mut self, names: Vec<&str>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows = self.row_suffix_to_indices(names)?;
        let rows: Vec<i32> = (0..self.records.len())
            .filter(|i| !rows.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();
        self.remove_records(rows)
    }

    fn concat(&mut self, others: Vec<&BaseAlignment>) -> PyResult<()> {
        check_empty_alignment(self)?;
        if let Err(_) = check_empty_list(&others) {
            return Ok(())
        };
        for aln in others.iter() {
            check_length_match(&self.records, &aln.records)?;
            check_length_match_i32(self.records[0].chunk_size,
                                   aln.records[0].chunk_size)?;
            for j in 0..self.records.len() {
                self.records[j].sequence.extend_from_slice(
                    &aln.records[j].sequence);
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
            "empty alignment"))
    }
    Ok(())
}

pub fn check_empty_list<T>(list: &Vec<T>) -> PyResult<()> {
    if list.len() == 0 {
        return Err(exceptions::ValueError::py_err(
            "empty list"))
    }
    Ok(())
}


pub fn check_row_index(aln_ptr: &BaseAlignment, i: usize) -> PyResult<()> {
    if aln_ptr.nrows()? <= i as i32 {
        return Err(exceptions::IndexError::py_err(
            "row index out of range"))
    }
    Ok(())
}

pub fn check_col_index(aln_ptr: &BaseAlignment, i: usize) -> PyResult<()> {
    if aln_ptr.ncols()? <= i as i32 {
        return Err(exceptions::IndexError::py_err(
            "column index out of range"))
    }
    Ok(())
}

pub fn check_chunk_size(aln_ptr: &BaseAlignment, list_ptr: &Vec<&str>)
-> PyResult<()> {
    if list_ptr.len() > 0 &&
        list_ptr[0].chars().count() != aln_ptr.chunk_size as usize {
        return Err(exceptions::ValueError::py_err(
            "chunk sizes do not match"))
    }
    Ok(())
}

pub fn check_length_match<T, U>(v1: &Vec<T>, v2: &Vec<U>) -> PyResult<()> {
    if v1.len() != v2.len() {
        return Err(exceptions::ValueError::py_err(
            "length mismatch"))
    }
    Ok(())
}

pub fn check_length_match_i32(len1: i32, len2:i32) -> PyResult<()> {
    if len1 != len2 {
        return Err(exceptions::ValueError::py_err(
            "length mismatch."))
    }
    Ok(())
}

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BaseAlignment>()?;

    Ok(())
}
