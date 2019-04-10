use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};

#[pyclass]
#[derive(Clone)]
pub struct SeqMatrix {
    pub data: Vec<String>,
    rows: usize,
    cols: usize,
}

pub fn new_seqmatrix(sequences: Vec<String>) -> Result<SeqMatrix, String> {
    let data = sequences.clone();
    let rows = data.len();    
    let cols = if rows > 0 { data[0].chars().count() } else { 0 };
    // Check whether each row has the same number of chars as the first row
    if rows > 0 {
        for row in sequences.iter() {
            if cols != row.chars().count() {
                return Err(format!("sequences have different character lengths"))
            }
        }
    }
    Ok(SeqMatrix{ data, rows, cols })
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
    pub fn _is_empty_matrix<'a>(&self) -> Result<(), String> {
        if self.rows == 0 {
            return Err("empty sequence matrix".to_owned())
        }
        Ok(())
    }

    /// Returns an error if a positive or negative index is greater than the size of the matrix
    pub fn _is_valid_row_index<'a>(&self, i: i32) -> Result<(), String> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.rows as i32 + i;
            if norm < 0 {
                return Err(format!("row ID ({}) is greater than the number of rows ({})", i, self.rows))
            }
        }
        // Check positive ID
        if i >= self.rows as i32 {
            return Err(format!("row ID ({}) is greater than the number of rows ({})", i, self.rows))
        }
        Ok(())
    }

    pub fn _is_valid_col_index<'a>(&self, i: i32) -> Result<(), String> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.cols as i32 + i;
            if norm < 0 {
                return Err(format!("column ID ({}) is greater than the number of column ({})", i, self.cols))
            }
        }
        // Check positive ID
        if i >= self.cols as i32 {
            return Err(format!("column ID is greater than the number of columns: {}", i))
        }
        Ok(())
    }

    // #endregion


    // Row methods
    // #region

    /// Returns a string sequence representing a row in the sequence matrix based on a given row index.
    pub fn _get_row<'a>(&self, id: i32) -> Result<String, String> {
        self._is_empty_matrix()?;
        // Convert negative index (count from end) to positive (count from start)
        let id: usize = if id < 0 { (self.rows as i32 + id) as usize } else { id as usize };
        Ok(self.data[id].to_string())
    }

    /// Returns a vector of string sequences representing rows in the sequence matrix based on the given vector of indices.
    pub fn _get_rows<'a>(&self, ids: Vec<i32>) -> Result<Vec<String>, String> {
        self._is_empty_matrix()?;
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        // Normalize row ids to positive values and get rows
        let result: Vec<String> = self._norm_rows(ids).into_iter()
            .map(|i| self.data[i].clone())
            .collect();
        Ok(result)
    }

    /// Removes rows from the sequence matrix based on a list of row indices.
    pub fn _remove_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_rows(ids, false)
    }

    /// Keep rows matching the specified row indices, and removes everything else.
    pub fn _retain_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_rows(ids, true)
    }

    /// Generalized method used to remove rows from the sequence matrix.
    fn _drop_rows<'a>(&mut self, ids: Vec<i32>, invert: bool) -> Result<(), String> {
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
        self.data = self.data.clone().into_iter().enumerate()
            .filter(|(i, _)| {
                if invert {
                    // rows in ids will be retained
                    rows.contains(i)
                } else {
                    // rows in ids will be removed
                    !rows.contains(i)
                }
                
            })
            .map(|(_, x)| x )
            .collect();
        Ok(())
    }

    /// Reorders rows based on a given ordered vector of row indices.
    pub fn _reorder_rows<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._is_empty_matrix()?;
        if ids.len() != self.rows {
            return Err(format!("number of ids ({}) is not equal to the number of rows ({})", ids.len(), self.rows))
        }
        if let Some(x) = ids.iter().max() {
            self._is_valid_row_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_row_index(*x)?;
        }
        
        // Normalize row ids to positive ids        
        let rows: Vec<usize> = self._norm_rows(ids);
        // Reorder using normalized row ids
        self.data = rows.into_iter().map(|i| self.data[i].clone()).collect();
        Ok(())
    }

    // #endregion


    // Column methods
    // #region

    /// Returns a single contiguous n-char column of the sequence matrix as vector of String for a given column index and chunk size.
    pub fn _get_chunk<'a>(&self, id: i32, chunk_size: usize) -> Result<Vec<String>, String> {
        self._is_empty_matrix()?;
        self._is_valid_col_index(id)?;
        let col: usize = if id < 0 { (self.rows as i32 + id) as usize } else { id as usize };
        let sequences: Vec<String> = self.data.iter()
            .map(|row| {
                let row: Vec<char> = row.chars().collect();
                let seq: String = row[col..col+chunk_size].iter().collect();
                seq
            })
            .collect();
        Ok(sequences)
    }

    /// Returns one or more contiguous n-char columns of the sequence matrix as vector of vector of String for a given vector of column indices and a chunk size.
    pub fn _get_chunks<'a>(&self, ids: Vec<i32>, chunk_size: usize) -> Result<Vec<Vec<String>>, String> {
        self._is_empty_matrix()?;
        let mut sorted_ids: Vec<i32> = ids.clone();
        sorted_ids.sort_unstable();
        if sorted_ids.len() == 0 {
            return Ok(vec![Vec::new()])
        } else if sorted_ids.len() == 1 {
            self._is_valid_col_index(sorted_ids[0])?;
        } else {
            self._is_valid_col_index(sorted_ids[0])?;
            self._is_valid_col_index(sorted_ids[sorted_ids.len()-1])?;
        }
        let seq_vec: Vec<Vec<char>> = self.data.iter()
            .map(|row| row.chars().collect())
            .collect();
        let sequences_vec: Vec<Vec<String>> = self._norm_cols(ids).into_iter()
            .map(|col| {
                let sequences: Vec<String> = seq_vec.iter()
                    .map(|row| row[col..col+chunk_size].iter().collect())
                    .collect();
                sequences
            })
            .collect();
        Ok(sequences_vec)
    }

    /// Returns a vector of string sequence representing a column in the sequence matrix based on the given index.
    pub fn _get_col<'a>(&self, id: i32) -> Result<Vec<String>, String> {
        self._get_chunk(id, 1)
    }

    /// Returns a vector of vector of string sequences representing columns in the sequence matrix based on the given vector of indices.
    pub fn _get_cols<'a>(&self, ids: Vec<i32>) -> Result<Vec<Vec<String>>, String> {
        self._get_chunks(ids, 1)
    }

    /// Removes columns from the sequence matrix based on a list of column indices.
    pub fn _remove_cols<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_cols(ids, false)
    }

    /// Keep columns matching the specified columns indices, and removes everything else.
    pub fn _retain_cols<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_cols(ids, true)
    }

    /// Generalized method used to remove columns from the sequence matrix.
    fn _drop_cols<'a>(&mut self, ids: Vec<i32>, invert: bool) -> Result<(), String> {
        self._is_empty_matrix()?;
        let mut sorted_ids: Vec<i32> = ids.clone();
        sorted_ids.sort_unstable();
        if sorted_ids.len() == 0 {
            return Ok(())
        } else if sorted_ids.len() == 1 {
            self._is_valid_col_index(sorted_ids[0])?;
        } else {
            self._is_valid_col_index(sorted_ids[0])?;
            self._is_valid_col_index(sorted_ids[sorted_ids.len()-1])?;
        }
        let cols: Vec<usize> = self._norm_cols(ids);
        self.data = self.data.clone().into_iter()
            .map(|row| {
                let sequence: String = row.char_indices()
                    .filter(|(i, _)| {
                        if invert {
                            // cols in ids will be retained
                            cols.contains(i)
                        } else {
                            // cols in ids will be removed
                            !cols.contains(i)
                        }
                    })
                    .map(|(_, x)| x )
                    .collect();
                sequence
            })
            .collect();
        Ok(())
    }

    /// Reorders rows based on a given ordered vector of row indices.
    pub fn _reorder_cols<'a>(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._is_empty_matrix()?;
        if ids.len() != self.cols {
            return Err(format!("number of ids ({}) is not equal to the number of columns ({})", ids.len(), self.rows))
        }
        if let Some(x) = ids.iter().max() {
            self._is_valid_col_index(*x)?;
        }
        if let Some(x) = ids.iter().min() {
            self._is_valid_col_index(*x)?;
        }
        
        // Normalize col ids to positive ids        
        let cols: Vec<usize> = self._norm_cols(ids);
        // Reorder using normalized col ids
        self.data = self.data.clone().into_iter()
            .map(|row| {
                let seq_vec: Vec<char> = row.chars().collect();
                let sequence: String = cols.iter()
                    .map(|j| seq_vec[*j])
                    .collect();
                sequence
            })
            .collect();
        Ok(())
    }

    // #endregion


    // SeqMatrix methods

    /// Concatenates sequence matrices across columns, preserving the number of rows.
    pub fn _concat<'a>(&mut self, others: Vec<&SeqMatrix>) -> Result<SeqMatrix, String> {
        self._is_empty_matrix()?;
        if others.len() == 0 {
            return Ok(self._copy())
        } else {
            let rows_true: usize = others.iter()
                .map(|m| if self.rows == m.rows { 1 } else { 0 })
                .sum();
            if others.len() != rows_true {
                return Err(format!("number of rows of other matrices is not equal to {}", self.rows))
            }
        }
        let mut sq = self._copy();
        for aln in others.iter() {
            for j in 0..self.data.len() {
                sq.data[j].push_str(&aln.data[j]);
            }
        }
        Ok(sq)
    }

    // TODO: implement clone()

    pub fn _copy(&self) -> SeqMatrix {
        SeqMatrix{
            data: self.data.clone(),
            rows: self.rows,
            cols: self.cols,
        }
    }

    // Utility methods
    // #region

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
    // #endregion
}

// Wrappers for pyo3
#[pymethods]
impl SeqMatrix {
    #[new]
    /// Creates a new SeqMatrix object from a list of sequences.
    fn __new__(obj: &PyRawObject, sequences: Vec<String>) -> PyResult<()> {
        let seq_matrix = match new_seqmatrix(sequences) {
            Ok(x) => x,
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        };
        // Instantiates the struct
        obj.init(|_| seq_matrix)
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
    // #region

    /// get_row(id, /)
    /// --
    /// 
    /// Returns a string sequence from the sequence matrix based on the given row index.
    fn get_row(&self, id: i32) -> PyResult<String> {
        match self._get_row(id) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// get_rows(ids, /)
    /// --
    /// 
    /// Returns a list of string sequences from the sequence matrix based on the given list of row indices.
    fn get_rows(&self, ids: Vec<i32>) -> PyResult<Vec<String>> {
        match self._get_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// remove_rows(ids, /)
    /// --
    /// 
    /// Removes rows from the sequence matrix based on a list of row indices.
    fn remove_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._remove_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// retain_records(indices, /)
    /// 
    /// Keep rows matching the specified row indices, and removes everything else.
    fn retain_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._retain_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// reorder_records(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    pub fn reorder_rows(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._reorder_rows(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }
    
    // #endregion


    // Column methods
    // #region

    /// get_chunk(id, chunk_size, /)
    /// --
    /// 
    /// Returns a single contiguous n-char column of the sequence matrix as list of str for a given column index and chunk size.
    fn get_chunk(&self, id: i32, chunk_size: usize) 
    -> PyResult<Vec<String>> {
        match self._get_chunk(id, chunk_size) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// get_chunks(ids, chunk_size, /)
    /// --
    /// 
    /// Returns one or more contiguous n-char columns of the sequence matrix as list of list of str for a given vector of column indices and a chunk size.
    fn get_chunks(&self, ids: Vec<i32>, chunk_size: usize) 
    -> PyResult<Vec<Vec<String>>> {
        match self._get_chunks(ids, chunk_size) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// get_col(id, /)
    /// --
    /// 
    /// Returns a list of sequence representing a column in the sequence matrix based on the given column index.
    fn get_col(&self, id: i32) -> PyResult<Vec<String>> {
        match self._get_col(id) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// get_cols(ids, /)
    /// --
    /// 
    /// Returns a list of list of sequences representing columns in the sequence matrix based on the given list of column indices.
    fn get_cols(&self, ids: Vec<i32>) -> PyResult<Vec<Vec<String>>> {
        match self._get_cols(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// remove_cols(indices, /)
    /// --
    /// 
    /// Removes many alignment columns simulatenously based on a list of column indices.
    pub fn remove_cols(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._remove_cols(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// retain_cols(indices, /)
    /// 
    /// Keep  alignment columns at the specified column indices and removes everything else.
    pub fn retain_cols(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._retain_cols(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    /// reorder_cols(ids, /)
    /// --
    /// 
    /// Reorders the alignment columns inplace based on a list of current column indices.
    pub fn reorder_cols(&mut self, ids: Vec<i32>) -> PyResult<()> {
        match self._reorder_cols(ids) {
            Ok(res) => Ok(res),
            Err(x) => Err(exceptions::IndexError::py_err(x)),
        }
    }

    // #endregion


    // SeqMatrix methods
    // #region

    /// concat(others, /)
    /// --
    /// 
    /// Returns a new sequence matrix by concatenating this matrix to a list of other matrices.
    pub fn concat(&mut self, others: Vec<&SeqMatrix>) -> PyResult<SeqMatrix> {
        match self._concat(others) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        }
    }

    /// copy()
    /// --
    /// 
    /// Returns a deep copy of the current sequence matrix.
    fn copy(&self) -> PyResult<SeqMatrix> {
        Ok(self._copy())
    }

    // #endregion


    // pub fn has(&self, query: &str, case_sensitive: bool, mode: &str, step_size: i32, chunk_size: i32) -> PyResult<Vec<i32>> {
    //     check_empty_alignment(self)?;
    //     let mut positions: Vec<i32> = Vec::new();
    //     let seq_vec: Vec<Vec<char>> = self.data.iter()
    //         .map(|seq| {seq.chars().collect()})
    //         .collect();
    //     let chunk_size = chunk_size as usize;
    //     let step_size = step_size as usize;
    //     if step_size < chunk_size {
    //         return Err(exceptions::ValueError::py_err(
    //                 "step_size is less than chunk_size"))
    //     }
    //     let query: String = match case_sensitive {
    //         true => query.to_string(),
    //         false => query.to_string().to_uppercase(),
    //     };
    //     for j in (0..seq_vec[0].len()).step_by(chunk_size) {
    //         let res: Vec<bool> = seq_vec.iter().map(
    //             |row| {
    //                 let mut chars: String = row[j..j+chunk_size]
    //                     .into_iter().collect();
    //                 if !case_sensitive {
    //                     chars = chars.to_uppercase();
    //                 }
    //                 if chars.contains(&query) {
    //                     true
    //                 } else {
    //                     false
    //                 }
    //             }
    //         ).collect();
    //         if mode == "any" {
    //             if res.contains(&true) {
    //                 let mut curr_positions: Vec<i32> = (j..j+chunk_size)
    //                     .map(|i| i as i32).collect();
    //                 positions.append(&mut curr_positions);
    //             }
    //         } else if mode == "all" {
    //             if !res.contains(&false) {
    //                 let mut curr_positions: Vec<i32> = (j..j+chunk_size)
    //                     .map(|i| i as i32).collect();
    //                 positions.append(&mut curr_positions);
    //             }
    //         } else {
    //             return Err(exceptions::ValueError::py_err(
    //                 "mode must be \"any\" or \"all\""))
    //         }
    //     }
    //     Ok(positions)
    // }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for SeqMatrix {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("SeqMatrix(nrows={nrows}, ncols={ncols})",
            nrows=self.nrows()?, ncols=self.ncols()?))
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(self.data.clone().join("\n"))
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

// #[pyfunction]
// pub fn from_list(sequences: Vec<String>) -> PyResult<BaseAlignment> {
//     let data: Vec<String> = sequences.iter().map(
//         |seq| seq.clone()
//     ).collect();
//     Ok(BaseAlignment{ data })
// }

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<SeqMatrix>()?;
    // m.add_function(wrap_function!(from_list))?;

    Ok(())
}
