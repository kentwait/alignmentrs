use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};
use std::fmt;

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
            let cnt = row.chars().count();
            if cols != cnt {
                return Err(format!("detected different sequences lengths: {} != {}", cols, cnt))
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
    pub fn _is_empty_matrix(&self) -> Result<(), String> {
        if self.rows == 0 {
            return Err("empty sequence matrix".to_owned())
        }
        Ok(())
    }

    /// Returns an error if a positive or negative index is greater than the size of the matrix
    pub fn _is_valid_row_index(&self, i: i32) -> Result<(), String> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.rows as i32 + i;
            if norm < 0 {
                return Err(format!("row ID ({}) is out of range [0,{})", i, self.rows))
            }
        }
        // Check positive ID
        if i >= self.rows as i32 {
            return Err(format!("row ID ({}) is out of range [0,{})", i, self.rows))
        }
        Ok(())
    }

    pub fn _is_valid_col_index(&self, i: i32) -> Result<(), String> {
        // Error if a negative index normalized to the start (0) is a negative number.
        if i < 0 {
            let norm = self.cols as i32 + i;
            if norm < 0 {
                return Err(format!("column ID ({}) is out of range [0,{})", i, self.cols))
            }
        }
        // Check positive ID
        if i >= self.cols as i32 {
            return Err(format!("column ID ({}) is out of range [0,{})", i, self.cols))
        }
        Ok(())
    }

    // #endregion


    // Row methods
    // #region

    /// Returns a string sequence representing a row in the sequence matrix based on a given row index.
    pub fn _get_row(&self, id: i32) -> Result<String, String> {
        self._is_empty_matrix()?;
        // Convert negative index (count from end) to positive (count from start)
        let id: usize = if id < 0 { (self.rows as i32 + id) as usize } else { id as usize };
        Ok(self.data[id].to_string())
    }

    /// Returns a vector of string sequences representing rows in the sequence matrix based on the given vector of indices.
    pub fn _get_rows(&self, ids: Vec<i32>) -> Result<Vec<String>, String> {
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
    pub fn _remove_rows(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_rows(ids, false)
    }

    /// Keep rows matching the specified row indices, and removes everything else.
    pub fn _retain_rows(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_rows(ids, true)
    }

    /// Generalized method used to remove rows from the sequence matrix.
    fn _drop_rows(&mut self, ids: Vec<i32>, invert: bool) -> Result<(), String> {
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
        self.rows = self.data.len();
        Ok(())
    }

    /// Reorders rows based on a given ordered vector of row indices.
    pub fn _reorder_rows(&mut self, ids: Vec<i32>) -> Result<(), String> {
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
    pub fn _get_chunk(&self, id: i32, chunk_size: usize) -> Result<Vec<String>, String> {
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
    pub fn _get_chunks(&self, ids: Vec<i32>, chunk_size: usize) -> Result<Vec<Vec<String>>, String> {
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
    pub fn _get_col(&self, id: i32) -> Result<Vec<String>, String> {
        self._get_chunk(id, 1)
    }

    /// Returns a vector of vector of string sequences representing columns in the sequence matrix based on the given vector of indices.
    pub fn _get_cols(&self, ids: Vec<i32>) -> Result<Vec<Vec<String>>, String> {
        self._get_chunks(ids, 1)
    }

    /// Removes columns from the sequence matrix based on a list of column indices.
    pub fn _remove_cols(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_cols(ids, false)
    }

    /// Keep columns matching the specified columns indices, and removes everything else.
    pub fn _retain_cols(&mut self, ids: Vec<i32>) -> Result<(), String> {
        self._drop_cols(ids, true)
    }

    /// Generalized method used to remove columns from the sequence matrix.
    fn _drop_cols(&mut self, ids: Vec<i32>, invert: bool) -> Result<(), String> {
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
        self.cols = self.data[0].len();
        Ok(())
    }

    /// Reorders rows based on a given ordered vector of row indices.
    pub fn _reorder_cols(&mut self, ids: Vec<i32>) -> Result<(), String> {
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
    pub fn _concat(&mut self, others: Vec<&SeqMatrix>) -> Result<(), String> {
        self._is_empty_matrix()?;
        if others.len() == 0 {
            return Ok(())
        } else {
            let rows_true: usize = others.iter()
                .map(|m| if self.rows == m.rows { 1 } else { 0 })
                .sum();
            if others.len() != rows_true {
                return Err(format!("number of rows of other matrices is not equal to {}", self.rows))
            }
        }
        // let mut sq = self._copy();
        for aln in others.iter() {
            for j in 0..self.data.len() {
                self.data[j].push_str(&aln.data[j]);
            }
        }
        self.cols = self.data[0].len();
        Ok(())
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
    pub fn _invert_rows(&self, ids: Vec<usize>) -> Vec<usize> {
        let rows: Vec<usize> = (0..self.rows)
                .filter(|i| !ids.contains(i) )
                .collect();
        rows
    }

    /// Returns column indices not found in the given vector of column indices.
    pub fn _invert_cols(&self, ids: Vec<usize>) -> Vec<usize> {
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
    pub fn concat(&mut self, others: Vec<&SeqMatrix>) -> PyResult<()> {
        match self._concat(others) {
            Ok(res) => Ok(res),
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        }
    }

    /// invert_rows(ids, /)
    /// --
    /// 
    /// Returns row indices that are not part of the given list of row indices.
    fn invert_rows(&self, ids: Vec<usize>) -> PyResult<Vec<usize>> {
        Ok(self._invert_rows(ids))
    }

    /// invert_cols(ids, /)
    /// --
    /// 
    /// Returns column indices that are not part of the given list of column indices.
    fn invert_cols(&self, ids: Vec<usize>) -> PyResult<Vec<usize>> {
        Ok(self._invert_cols(ids))
    }

    /// copy()
    /// --
    /// 
    /// Returns a deep copy of the current sequence matrix.
    fn copy(&self) -> PyResult<SeqMatrix> {
        Ok(self._copy())
    }
    // #endregion
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

// Implements equality comparison between SeqMatrix structs
impl PartialEq for SeqMatrix {
    fn eq(&self, other: &SeqMatrix) -> bool {
        self.data == other.data && self.rows == other.rows && self.cols == other.cols
    }
}

// Implements Debug in order to use format! and other printout methods
impl fmt::Debug for SeqMatrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SeqMatrix {{ data: {:?}, rows: {}, cols: {} }}", self.data, self.rows, self.cols)
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

mod test {
    use super::*;

    // Test SeqMatrix creation
    #[test]
    fn test_new_seqmatrix() {
        let exp = SeqMatrix{ 
            data: vec![
                "atcg".to_string(),
                "atgg".to_string(),
                "atcc".to_string(),
                "tagc".to_string(),
            ],
            rows: 4,
            cols: 4,
        };
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        assert_eq!(exp, res);
    }

    // Test SeqMatrix methods

    // Test getter methods
    #[test]
    fn test_nrows() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        assert_eq!(res.rows, 4);
    }

    #[test]
    fn test_ncols() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        assert_eq!(res.cols, 4);
    }

    // TODO: Place subsequent tests in another file

    // Test methods that used for checking
    #[test]
    fn test_is_empty_matrix() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        res._is_empty_matrix().unwrap();
    }

    #[test]
    #[should_panic(expected = "empty sequence matrix")]
    fn test_is_empty_matrix_empty() {
        let res = new_seqmatrix(vec![]).unwrap();
        
        res._is_empty_matrix().unwrap();
    }

    #[test]
    fn test_is_valid_row_index() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        res._is_valid_row_index(3).unwrap();
    }

    #[test]
    #[should_panic(expected = "row ID (4) is out of range [0,4)")]
    fn test_is_valid_row_index_invalid() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        res._is_valid_row_index(4).unwrap();
    }

    #[test]
    fn test_is_valid_col_index() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        res._is_valid_col_index(3).unwrap();
    }

    #[test]
    #[should_panic(expected = "column ID (4) is out of range [0,4)")]
    fn test_is_valid_col_index_invalid() {
        let res = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        res._is_valid_col_index(4).unwrap();
    }

    // Test row methods
    // Test _get_row and _get_rows
    #[test]
    fn test_get_row() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_row(0).unwrap();
        assert_eq!(res, "atcg");
    }

    #[test]
    fn test_get_row_negative_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_row(-1).unwrap();
        assert_eq!(res, "tagc");
    }

    #[test]
    fn test_get_rows() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_rows(vec![0,2]).unwrap();
        assert_eq!(res, vec!["atcg", "atcc"]);
    }

    #[test]
    fn test_get_rows_mixed_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_rows(vec![0,-1]).unwrap();
        assert_eq!(res, vec!["atcg", "tagc"]);
    }

    #[test]
    // Test remove rows
    fn test_remove_rows() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec!["atgg".to_string(), "tagc".to_string()]).unwrap();

        mat._remove_rows(vec![0, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    #[test]
    // Test retain rows
    fn test_retain_rows() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec!["atcg".to_string(), "atcc".to_string()]).unwrap();

        mat._retain_rows(vec![0, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    #[test]
    // Test retain rows base - drop rows, true
    fn test_drop_rows_true() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec!["atcg".to_string(), "atcc".to_string()]).unwrap();

        mat._drop_rows(vec![0, 2], true).unwrap();
        assert_eq!(mat, exp);
    }

    #[test]
    // Test remove rows base - drop rows, false
    fn test_drop_rows_false() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec!["atgg".to_string(), "tagc".to_string()]).unwrap();

        mat._drop_rows(vec![0, 2], false).unwrap();
        assert_eq!(mat, exp);
    }

    #[test]
    // Test reorder rows
    fn test_reorder_rows() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "atgg".to_string(),  // 1
            "tagc".to_string(),  // 3
            "atcg".to_string(),  // 0
            "atcc".to_string(),  // 2
        ]).unwrap();

        mat._reorder_rows(vec![1, 3, 0, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    // Test column methods
    // Test _get_chunk and _get_chunks
    #[test]
    fn test_get_chunk() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunk(0, 1).unwrap();
        assert_eq!(res, vec!["a","a","a","t"]);
    }

    #[test]
    fn test_get_chunk_negative_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunk(-1, 1).unwrap();
        assert_eq!(res, vec!["g","g","c","c"]);
    }

    #[test]
    fn test_get_chunk_3() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunk(0, 3).unwrap();
        assert_eq!(res, vec!["atc","atg","atc","tag"]);
    }

    #[test]
    fn test_get_chunks() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunks(vec![0,2], 1).unwrap();
        assert_eq!(res, vec![vec!["a","a","a","t"],vec!["c","g","c","g"]]);
    }

    #[test]
    fn test_get_chunks_mixed_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunks(vec![0,-1], 1).unwrap();
        assert_eq!(res, vec![vec!["a","a","a","t"],vec!["g","g","c","c"]]);
    }

    #[test]
    fn test_get_chunks_3() {
        let mat = new_seqmatrix(vec![
            "atcgt".to_string(),
            "atggt".to_string(),
            "atccg".to_string(),
            "tagcc".to_string(),
        ]).unwrap();
        
        let res = mat._get_chunks(vec![0,2], 3).unwrap();
        assert_eq!(res, vec![vec!["atc","atg","atc","tag"],vec!["cgt","ggt","ccg","gcc"]]);
    }

    // Tests _get_col and _get_cols
    #[test]
    fn test_get_col() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_col(0).unwrap();
        assert_eq!(res, vec!["a","a","a","t"]);
    }

    #[test]
    fn test_get_col_negative_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_col(-1).unwrap();
        assert_eq!(res, vec!["g","g","c","c"]);
    }

    #[test]
    fn test_get_cols() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_cols(vec![0,2]).unwrap();
        assert_eq!(res, vec![vec!["a","a","a","t"],vec!["c","g","c","g"]]);
    }

    #[test]
    fn test_get_cols_mixed_index() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        
        let res = mat._get_cols(vec![0,-1]).unwrap();
        assert_eq!(res, vec![vec!["a","a","a","t"],vec!["g","g","c","c"]]);
    }

    // Tests _remove_cols
    #[test]
    fn test_remove_cols() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "tg".to_string(),
            "tg".to_string(),
            "tc".to_string(),
            "ac".to_string(),
        ]).unwrap();

        mat._remove_cols(vec![0, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    // Tests _retain_cols
    #[test]
    fn test_retain_cols() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "ac".to_string(),
            "ag".to_string(),
            "ac".to_string(),
            "tg".to_string(),
        ]).unwrap();

        mat._retain_cols(vec![0, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    // Tests _retain_cols base - _drop_cols, true
    #[test]
    fn test_drop_cols_true() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "ac".to_string(),
            "ag".to_string(),
            "ac".to_string(),
            "tg".to_string(),
        ]).unwrap();

        mat._drop_cols(vec![0, 2], true).unwrap();
        assert_eq!(mat, exp);
    }

    // Test _remove_cols base - _drop_cols, false
    #[test]
    fn test_drop_cols_false() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "tg".to_string(),
            "tg".to_string(),
            "tc".to_string(),
            "ac".to_string(),
        ]).unwrap();

        mat._drop_cols(vec![0, 2], false).unwrap();
        assert_eq!(mat, exp);
    }

    // Test _reorder_cols
    #[test]
    fn test_reorder_cols() {
        let mut mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "gatc".to_string(),
            "gatg".to_string(),
            "catc".to_string(),
            "ctag".to_string(),
        ]).unwrap();

        mat._reorder_cols(vec![3, 0, 1, 2]).unwrap();
        assert_eq!(mat, exp);
    }

    // Test normalization of index values
    #[test]
    fn test_norm_rows() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        assert_eq!(mat._norm_rows(vec![0, -1, 2, -3, -4, 3]), vec![0, 3, 2, 1, 0, 3]);
    }

    #[test]
    fn test_norm_cols() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        assert_eq!(mat._norm_cols(vec![0, -1, 2, -3, -4, 3]), vec![0, 3, 2, 1, 0, 3]);
    }

    // Test index inversion
    #[test]
    fn test_invert_rows() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        assert_eq!(mat._invert_rows(vec![3, 0]), vec![1, 2]);
    }

    #[test]
    fn test_invert_cols() {
        let mat = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();
        assert_eq!(mat._invert_cols(vec![3, 0]), vec![1, 2]);
    }

    // Tests concatenation
    #[test]
    fn test_concat() {
        let mut mat1 = new_seqmatrix(vec![
            "atcg".to_string(),
            "atgg".to_string(),
            "atcc".to_string(),
            "tagc".to_string(),
        ]).unwrap();

        let mat2 = new_seqmatrix(vec![
            "aaaa".to_string(),
            "tttt".to_string(),
            "cccc".to_string(),
            "gggg".to_string(),
        ]).unwrap();

        let exp = new_seqmatrix(vec![
            "atcgaaaa".to_string(),
            "atggtttt".to_string(),
            "atcccccc".to_string(),
            "tagcgggg".to_string(),
        ]).unwrap();

        mat1._concat(vec![&mat2]).unwrap();
        assert_eq!(mat1, exp);
    }
}