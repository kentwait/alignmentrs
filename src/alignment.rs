use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};

use crate::record::Record;


#[pyclass]
#[derive(Clone)]
pub struct SeqMatrix {
    pub data: Vec<String>,
    pub rows: usize,
    pub cols: usize,
}

// Rust functions
impl SeqMatrix {

    // Getter functions
    // #region

    /// Returns the number of rows in the sequence matrix.
    pub fn _nrows(&self) -> usize {
        self.rows
    }

    /// Returns the number of columns in the sequence matrix.
    pub fn _ncols(&self) -> usize {
        self.cols
    }

    // #endregion


    // Error methods
    // #region

    /// Returns an error if the matrix is empty.
    pub fn _is_empty_matrix<'a>(&self) -> Result<(), &'a str> {
        if self.rows == 0 {
            return Err("empty sequence matrix")
        }
        Ok(())
    }

    /// Returns an error if a positive or negative index is greater than the size of the matrix
    pub fn _is_valid_row_index<'a>(&self, i: i32) -> Result<(), &'a str> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.rows as i32 + i;
            if norm < 0 {
                return Err(&format!("row ID ({}) is greater than the number of rows ({})", i, self.rows))
            }
        }
        // Check positive ID
        if i as usize >= self.rows {
            return Err(&format!("row ID ({}) is greater than the number of rows ({})", i, self.rows))
        }
        Ok(())
    }

    pub fn _is_valid_col_index<'a>(&self, i: i32) -> Result<(), &'a str> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.cols as i32 + i;
            if norm < 0 {
                return Err(&format!("column ID ({}) is greater than the number of column ({})", i, self.cols))
            }
        }
        // Check positive ID
        if i as usize >= self.cols {
            return Err(&format!("column ID is greater than the number of columns: {}", i))
        }
        Ok(())
    }

    // #endregion


    // Row methods
    // #region

    /// Returns a string sequence from the sequence matrix based on the given index.
    pub fn _get_row<'a>(&self, id: i32) -> Result<String, &'a str> {
        self._is_empty_matrix()?;
        // Convert negative index (count from end) to positive (count from start)
        let id: usize = if id < 0 { (self.rows as i32 + id) as usize } else { id as usize };
        Ok(self.data[id].to_string())
    }

    /// Returns a vector of string sequences from the sequence matrix based on the given vector of indices.
    pub fn _get_rows<'a>(&self, ids: Vec<i32>) -> Result<Vec<String>, &'a str> {
        self._is_empty_matrix()?;
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        // Normalize row ids to positive values and get rows
        let result: Vec<String> = self._norm_rows(ids).into_iter()
            .map(|i| self.data[i])
            .collect();
        Ok(result)
    }

    /// Removes rows from the sequence matrix based on a list of row indices.
    pub fn _remove_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), &'a str> {
        self._is_empty_matrix()?;
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        
        // Normalize row ids to positive ids
        let rows: Vec<usize> = self._norm_rows(ids);
        // Keep data whose index is not found in the rows vector
        // Remove if index is in the rows vector
        self.data = self.data.into_iter().enumerate()
            .filter(|(i, _)| !rows.contains(i))
            .map(|(_, x)| x )
            .collect();
        Ok(())
    }

    /// Keep rows matching the specified row indices, and removes everything else.
    pub fn _retain_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), &'a str> {
        self._is_empty_matrix()?;
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        
        // Normalize row ids to positive ids        
        let rows: Vec<usize> = self._norm_rows(ids);
        // Keep data whose index is in the rows vector
        // Remove if index is not found in the rows vector
        self.data = self.data.into_iter().enumerate()
            .filter(|(i, _)| rows.contains(i))
            .map(|(_, x)| x )
            .collect();
        Ok(())
    }

    pub fn _reorder_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), &'a str> {
        self._is_empty_matrix()?;
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        
        // Normalize row ids to positive ids        
        let rows: Vec<usize> = self._norm_rows(ids);
        // Reorder using normalized row ids
        self.data = rows.into_iter().map(|i| self.data[i]).collect();
        Ok(())
    }

    // #endregion


    // Column methods




    // Utility methods

    /// Converts row indices into positive-value row indices.
    pub fn _norm_rows(&self, ids: Vec<i32>) -> Vec<usize> {
        let normed_rows: Vec<usize> = ids.iter()
            .map(|i| {
                if *i >= 0 { 
                    *i as usize
                } else {
                    (self.rows as i32 + *i) as usize
                }
            })
            .collect();
        normed_rows
    }

    /// Converts column indices into positive-value column indices.
    pub fn _norm_cols(&self, ids: Vec<i32>) -> Vec<usize> {
        let normed_cols: Vec<usize> = ids.iter()
            .map(|i| {
                if *i >= 0 { 
                    *i as usize
                } else {
                    (self.cols as i32 + *i) as usize
                }
            })
            .collect();
        normed_cols
    }

    /// Returns row indices not found in the given vector of row indices.
    pub fn _invert_rows<'a>(&self, ids: Vec<usize>) -> Vec<usize> {
        let rows: Vec<usize> = (0..self.rows)
                .filter(|i| !ids.contains(i) )
                .collect();
        rows
    }

    /// Returns column indices not found in the given vector of column indices.
    pub fn _invert_cols<'a>(&self, ids: Vec<usize>) -> Vec<usize> {
        let cols: Vec<usize> = (0..self.cols)
                .filter(|i| !ids.contains(i) )
                .collect();
        cols
    }
}

// Wrappers for pyo3
#[pymethods]
impl SeqMatrix {
    #[new]
    /// Creates a new SeqMatrix object from a list of sequences.
    fn __new__(obj: &PyRawObject, sequences: Vec<String>) -> PyResult<()> {
        let rows: usize = sequences.len();
        let mut cols: usize = 0;
        // Check whether each row has the same number of chars as the first row
        if sequences.len() > 0 {
            cols = sequences[0].chars().count();
            for row in sequences.iter() {
                if cols != row.chars().count() {
                    return Err(exceptions::ValueError::py_err(
                    "sequences have different character lengths"))
                }
            }
        }
        let data = sequences.clone();
        // Instantiates the struct
        obj.init(|_| {
            SeqMatrix{ data, rows, cols }
        })
    }
    
    #[getter]
    /// int: Returns the number of rows in the BaseAlignment.
    fn nrows(&self) -> PyResult<i32> {
        Ok(self._nrows() as i32)
    }

    #[getter]
    /// int: Returns the number of columns in the alignment.
    fn ncols(&self) -> PyResult<i32> {
        Ok(self._ncols() as i32)
    }

    #[getter]
    /// list of str: Returns the list of sequences.
    fn sequences(&self) -> PyResult<Vec<String>> {
        Ok(self.data.clone())
    }

    // Row methods

    /// get_row(id, /)
    /// --
    /// 
    /// Returns a string sequence from the sequence matrix based on the given row index.
    fn get_row(&self, id: i32) -> PyResult<String> {
        match self._get_row(id) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// get_rows(ids, /)
    /// --
    /// 
    /// Returns a list of string sequences from the sequence matrix based on the given list of row indices.
    fn get_rows(&self, ids: Vec<i32>) -> PyResult<Vec<String>> {
        match self._get_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// remove_rows(ids, /)
    /// --
    /// 
    /// Removes rows from the sequence matrix based on a list of row indices.
    fn remove_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._remove_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// retain_records(indices, /)
    /// 
    /// Keep rows matching the specified row indices, and removes everything else.
    fn retain_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._retain_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// reorder_records(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    pub fn reorder_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._reorder_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// insert_records(position, ids, descriptions, sequences, /)
    /// 
    /// Inserts one or more samples at the specified position.
    // fn insert_records(&mut self, mut row: i32, values: Vec<&Record>) -> PyResult<()> {
    //     if self.nrows()? < row {
    //         return Err(exceptions::IndexError::py_err(
    //             "row index out of range."))
    //     }
    //     for value in values.into_iter() {
    //         self.records.insert(row as usize, value.clone());
    //         row += 1;
    //     }
    //     Ok(())
    // }    

    // fn replace_record(&mut self, row: i32, value: &Record) -> PyResult<()> {
    //     self.replace_records(vec![row], vec![value])
    // }

    // fn replace_records(&mut self, rows: Vec<i32>, values: Vec<&Record>) -> PyResult<()> {
    //     check_length_match(&rows, &values)?;
    //     check_empty_alignment(self)?;
    //     if let Some(x) = rows.iter().max() {
    //         check_row_index(self, *x as usize)?;
    //     }
    //     for (i, row) in rows.into_iter().map(|x| x as usize).enumerate() {
    //         // TODO: Make function that checks if records have the same length and string length
    //         check_length_match_i32(values[i].len_str()?, self.records[i].len_str()?)?;
    //         check_length_match_i32(values[i].len()?, self.records[i].len()?)?;
    //         self.records[row] = values[i].clone();
    //     }
    //     Ok(())
    // }

    // pub fn drain_rows(&mut self, mut rows: Vec<i32>) -> PyResult<BaseAlignment> {
    //     check_empty_alignment(self)?;
    //     rows.sort_unstable();
    //     rows.dedup();
    //     if let Some(x) = rows.iter().max() {
    //         check_row_index(self, *x as usize)?;
    //     }
    //     let mut data: Vec<String> = Vec::new();
    //     let mut i: i32 = 0;
    //     while i != self.data.len() as i32 {
    //         if rows.contains(&i) {
    //             data.push(self.data.remove(i as usize));
    //         } else {
    //             i += 1;
    //         }
    //     }
    //     Ok(BaseAlignment{ data }) 
    // }
    



    // Column methods

    pub fn get_col(&self, col: i32) -> PyResult<Vec<String>> {
        self.get_chunk(col, 1)
    }

    /// get_cols(col_indices, /)
    /// --
    /// 
    /// Returns a new RawAlignment object containing the specific
    /// alignment columns based on a list of indices. 
    pub fn get_cols(&self, cols: Vec<i32>) -> PyResult<Vec<Vec<String>>> {
        self.get_chunks(cols, 1)
    }

    pub fn get_chunk(&self, col: i32, chunk_size: i32) 
    -> PyResult<Vec<String>> {
        check_empty_alignment(self)?;
        check_col_index(self, col as usize)?;
        let col = col as usize;
        let chunk_size = chunk_size as usize;
        let seq_vec: Vec<Vec<char>> = self.data.iter()
            .map(|seq| {seq.chars().collect()})
            .collect();
        let sequences: Vec<String> = seq_vec.into_iter()
            .map(|row| row[col..col+chunk_size].into_iter().collect())
            .collect();
        Ok(sequences)
    }

    pub fn get_chunks(&self, cols: Vec<i32>, chunk_size: i32) 
    -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        let chunk_size = chunk_size as usize;
        let seq_vec: Vec<Vec<char>> = self.data.iter()
            .map(|seq| seq.chars().collect())
            .collect();
        let sequences_vec: Vec<Vec<String>> = cols.into_iter()
            .map(|col| {
                let col = col as usize;
                let sequences: Vec<String> = seq_vec.iter()
                    .map(|row| row[col..col+chunk_size].iter().collect())
                    .collect();
                sequences
            })
            .collect();
        Ok(sequences_vec)
    }

    // insert
    // fn insert_col(&mut self, col: usize, value: Vec<&str>) -> PyResult<()> {
    //     self.insert_cols(col, vec![value])
    // }

    // fn insert_cols(&mut self, mut col: usize, values: Vec<Vec<&str>>)
    // -> PyResult<()> {
    //     if values.len() == 0 {
    //         return Ok(())
    //     }
    //     // TODO: Add immediate return to other methods when value length is 0
    //     check_col_index(self, col)?;
    //     check_length_match(&values[0], &self.data)?;
    //     check_chunk_size(self, &values[0])?;
    //     for value in values.iter() {
    //         for i in 0..self.records.len() {
    //             self.records[i].sequence.insert(
    //                 col as usize, value[i].to_string());
    //         }
    //         col += 1;
    //     }
    //     Ok(())
    // }

    // for value in values.into_iter() {
    //     self.records.insert(row as usize, value.clone());
    //     row += 1;
    // }
    // Ok(())

    // append
    // fn append_col(&mut self, value: Vec<&str>) -> PyResult<()> {
    //     self.append_cols(vec![value])
    // }

    // fn append_cols(&mut self, values: Vec<Vec<&str>>) -> PyResult<()> {
    //     if values.len() == 0 {
    //         return Ok(())
    //     }
    //     check_length_match(&values[0], &self.records)?;
    //     check_chunk_size(self, &values[0])?;
    //     for i in 0..self.records.len() {
    //         for value in values.iter() {
    //             self.records[i].sequence.push(value[i].to_string());
    //         }
    //     }
    //     Ok(())
    // }

    /// remove_cols(indices, /)
    /// --
    /// 
    /// Removes many alignment columns simulatenously based on a
    /// list of column indices.
    pub fn remove_cols(&mut self, mut cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        cols.sort_unstable();
        cols.dedup();
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        for row in 0..self.data.len() {
            let sequence: String = self.data[row].char_indices()
                .filter(|(i, _)| !cols.contains(&(*i as i32)))
                .map(|(_, x)| x )
                .collect();
            self.data[row] = sequence;
        }
        Ok(())
    }

    /// retain_cols(indices, /)
    /// 
    /// Keep  alignment columns at the specified column indices and
    /// removes everything else.
    pub fn retain_cols(&mut self, cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let cols: Vec<i32> = (0..self.ncols()?)
            .filter(|i| !cols.contains(&(*i as i32)) )
            .map(|i| i as i32 )
            .collect();        
        self.remove_cols(cols)
    }

    pub fn drain_cols(&mut self, cols: Vec<i32>) -> PyResult<BaseAlignment> {
        let mut aln = self.copy()?;
        aln.retain_cols(cols.clone())?;
        self.remove_cols(cols.clone())?;
        Ok(aln)
    }

    /// replace_col(coordinate, sequence, /)
    /// --
    /// 
    /// Replaces the all the characters in an alignment column.
    // fn replace_col(&mut self, col: i32, value: Vec<&str>) -> PyResult<()> {
    //     self.replace_cols(vec![col], vec![value])
    // }

    /// replace_cols(coordinates, sequences, /)
    /// --
    /// 
    /// Replaces the all the characters in each specified column from a list of
    /// alignment columns.
    // fn replace_cols(&mut self, cols: Vec<i32>, values: Vec<Vec<&str>>)
    // -> PyResult<()> {
    //     check_length_match(&cols, &values)?;
    //     check_empty_alignment(self)?;
    //     if let Some(x) = cols.iter().max() {
    //         check_col_index(self, *x as usize)?;
    //     }
    //     for (i, col) in cols.iter().enumerate() {
    //         check_length_match(&self.records, &values[i])?;
    //         for row in 0..self.records.len() {
    //             self.records[row].sequence[*col as usize] = values[i][row].to_string();
    //         }
    //     }
    //     Ok(())
    // }

    /// reorder_cols(ids, /)
    /// --
    /// 
    /// Reorders the alignment columns inplace based on a list of current
    /// column indices.
    pub fn reorder_cols(&mut self, cols: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        if cols.len() != self.data.len() {
            return Err(exceptions::ValueError::py_err(
                "length mismatch."))
        }
        for i in 0..self.data.len() {
            let seq_vec: Vec<char> = self.data[i].chars().collect();
            let sequence: String = cols.iter()
                .map(|j| seq_vec[(*j) as usize])
                .collect();
            self.data[i] = sequence;
        }
        Ok(())
    }

    pub fn has(&self, query: &str, case_sensitive: bool, mode: &str, step_size: i32, chunk_size: i32) -> PyResult<Vec<i32>> {
        check_empty_alignment(self)?;
        let mut positions: Vec<i32> = Vec::new();
        let seq_vec: Vec<Vec<char>> = self.data.iter()
            .map(|seq| {seq.chars().collect()})
            .collect();
        let chunk_size = chunk_size as usize;
        let step_size = step_size as usize;
        if step_size < chunk_size {
            return Err(exceptions::ValueError::py_err(
                    "step_size is less than chunk_size"))
        }
        let query: String = match case_sensitive {
            true => query.to_string(),
            false => query.to_string().to_uppercase(),
        };
        for j in (0..seq_vec[0].len()).step_by(chunk_size) {
            let res: Vec<bool> = seq_vec.iter().map(
                |row| {
                    let mut chars: String = row[j..j+chunk_size]
                        .into_iter().collect();
                    if !case_sensitive {
                        chars = chars.to_uppercase();
                    }
                    if chars.contains(&query) {
                        true
                    } else {
                        false
                    }
                }
            ).collect();
            if mode == "any" {
                if res.contains(&true) {
                    let mut curr_positions: Vec<i32> = (j..j+chunk_size)
                        .map(|i| i as i32).collect();
                    positions.append(&mut curr_positions);
                }
            } else if mode == "all" {
                if !res.contains(&false) {
                    let mut curr_positions: Vec<i32> = (j..j+chunk_size)
                        .map(|i| i as i32).collect();
                    positions.append(&mut curr_positions);
                }
            } else {
                return Err(exceptions::ValueError::py_err(
                    "mode must be \"any\" or \"all\""))
            }
        }
        Ok(positions)
    }

    /// subset(row_indices, column_indices, /)
    /// --
    /// 
    /// Returns the subset of rows and columns in the alignment as a new
    /// RawAlignment.
    // fn subset(&self, rows: Vec<i32>, cols: Vec<i32>)
    // -> PyResult<BaseAlignment> {
    //     check_empty_alignment(self)?;
    //     if let Some(x) = rows.iter().max() {
    //         check_row_index(self, *x as usize)?;
    //     }
    //     if let Some(x) = cols.iter().max() {
    //         check_col_index(self, *x as usize)?;
    //     }
    //     let mut records: Vec<Record> = Vec::new();
    //     let chunk_size = self.chunk_size;
    //     let ncols = cols.len();
    //     for row in rows.into_iter().map(|x| x as usize) {
    //         let sequence: Vec<String> = self.records[row].sequence.iter().enumerate()
    //             .filter(|(i, _)| cols.contains(&(*i as i32)) )
    //             .map(|(_, item)| item.to_string() )
    //             .collect();
    //         // Make sure sequence length == ncols
    //         if sequence.len() != ncols {
    //             return Err(exceptions::ValueError::py_err(
    //                 format!("Unexpected number of columns: {} != {}",
    //                     sequence.len(), ncols)))
    //         }
    //         records.push(Record{
    //             id: self.records[row].id.to_string(),
    //             description: self.records[row].description.to_string(),
    //             sequence,
    //             chunk_size,
    //         })
    //     }
    //     Ok(BaseAlignment{ records, chunk_size })
    // }


    pub fn concat(&mut self, others: Vec<&BaseAlignment>) -> PyResult<()> {
        check_empty_alignment(self)?;
        if let Err(_) = check_empty_list(&others) {
            return Ok(())
        };
        for aln in others.iter() {
            check_length_match(&self.data, &aln.data)?;
            for j in 0..self.data.len() {
                self.data[j].push_str(&aln.data[j]);
            }
        }
        Ok(())
    }

    pub fn copy(&self) -> PyResult<BaseAlignment> {
        Ok(BaseAlignment{ data: self.data.clone() })
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for SeqMatrix {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("BaseAlignment(nrows={nrows}, ncols={ncols})",
            nrows=self.nrows()?, ncols=self.ncols()?))
    }

    fn __str__(&self) -> PyResult<String> {
        if self.nrows()? == 0 {
            return Ok(String::new())
        }
        let mut fasta_strings: Vec<String> = Vec::new();
        for seq in self.data.iter() {
            fasta_strings.push(seq.to_string());
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

#[pyfunction]
pub fn from_list(sequences: Vec<String>) -> PyResult<BaseAlignment> {
    let data: Vec<String> = sequences.iter().map(
        |seq| seq.clone()
    ).collect();
    Ok(BaseAlignment{ data })
}

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
    m.add_class::<SeqMatrix>()?;
    m.add_function(wrap_function!(from_list))?;

    Ok(())
}
