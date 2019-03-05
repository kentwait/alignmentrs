use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use pyo3::class::gc::{PyGCProtocol, PyVisit, PyTraverseError};

use crate::record::Record;

#[pyclass]
#[derive(Clone)]
pub struct BaseAlignment {

    pub data: Vec<String>,

}

#[pymethods]
impl BaseAlignment {
    #[new]
    /// Creates a new BaseAlignment object from a list of ids, descriptions,
    /// and sequences.
    fn __new__(obj: &PyRawObject, records: Vec<&Record>) -> PyResult<()> {
        // TODO: Check if all vectors have the same size.
        let mut data: Vec<String> = Vec::with_capacity(records.len());
        for record in records.into_iter() {
            data.push(record.sequence.to_string());
        }
        obj.init(|_| {
            BaseAlignment{ data }
        })
    }
    
    #[getter]
    /// int: Returns the number of rows in the BaseAlignment.
    pub fn nrows(&self) -> PyResult<i32> {
        Ok(self.data.len() as i32)
    }

    #[getter]
    /// int: Returns the number of columns in the alignment.
    pub fn ncols(&self) -> PyResult<i32> {
        if self.data.len() == 0 {
            return Ok(0)
        }
        Ok(self.data[0].len() as i32)
    }

    #[getter]
    /// int: Returns the number of aligned characters (ncols * chunk_size).
    pub fn nchars(&self) -> PyResult<i32> {
        if self.data.len() == 0 {
            return Ok(0)
        }
        Ok(self.data[0].len() as i32)
    }

    #[getter]
    /// list of str: Returns the list of sequences.
    pub fn sequences(&self) -> PyResult<Vec<String>> {
        Ok(self.data.clone())
    }

    // Row methods
    pub fn get_row(&self, row: i32) -> PyResult<String> {
        check_empty_alignment(self)?;
        check_row_index(self, row as usize)?;
        Ok(self.data[row as usize].to_string())
    }


    pub fn get_rows(&self, rows: Vec<i32>) -> PyResult<Vec<String>> {
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut sequences: Vec<String> = Vec::new();
        for row in rows.into_iter().map(|x| x as usize) {
            sequences.push(self.data[row].to_string());
        }
        Ok(sequences)
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

    /// remove_records(indices, /)
    /// --
    /// 
    /// Removes many entries simulatenously based on a
    /// list of row indices.
    fn remove_rows(&mut self, mut rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        rows.sort_unstable();
        rows.dedup();
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        rows.reverse();
        for row in rows.iter().map(|x| *x as usize) {
            self.data.remove(row);
        }
        Ok(())
    }

    /// retain_records(indices, /)
    /// 
    /// Keep entries at the specified row indices, and removes
    /// everything else.
    fn retain_rows(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_empty_alignment(self)?;
        let rows: Vec<i32> = self.invert_rows(rows)?;
        self.remove_rows(rows)
    }

    fn drain_rows(&mut self, mut rows: Vec<i32>) -> PyResult<BaseAlignment> {
        check_empty_alignment(self)?;
        rows.sort_unstable();
        rows.dedup();
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut data: Vec<String> = Vec::new();
        let mut i: i32 = 0;
        while i != self.data.len() as i32 {
            if rows.contains(&i) {
                data.push(self.data.remove(i as usize));
            } else {
                i += 1;
            }
        }
        Ok(BaseAlignment{ data }) 
    }

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

    /// reorder_records(ids, /)
    /// --
    /// 
    /// Reorders the sequences inplace based on a list of current row indices.
    fn reorder_rows(&mut self, rows: Vec<i32>) -> PyResult<()> {
        check_length_match(&rows, &self.data)?;
        check_empty_alignment(self)?;
        if let Some(x) = rows.iter().max() {
            check_row_index(self, *x as usize)?;
        }
        let mut data: Vec<String> = Vec::with_capacity(rows.len());
        for row in rows.into_iter().map(|x| x as usize) {
            data.push(self.data[row].clone());
        }
        self.data = data;
        Ok(())
    }



    // Column methods

    pub fn get_col(&self, col: i32) -> PyResult<Vec<String>> {
        check_empty_alignment(self)?;
        check_col_index(self, col as usize)?;
        let mut sequences: Vec<String> = Vec::new();
        for seq in self.data.iter() {
            let seq_vec: Vec<char> = seq.chars().collect();
            sequences.push(seq_vec[col as usize].to_string());
        }
        Ok(sequences)
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
        let mut out: Vec<Vec<String>> = Vec::new();
        for i in cols.into_iter().map(|x| x as usize) {
            let mut sequences: Vec<String> = Vec::new();
            for seq in self.data.iter() {
                let seq_vec: Vec<char> = seq.chars().collect();
                sequences.push(seq_vec[i].to_string());
            }
            out.push(sequences);
        }
        Ok(out)
    }

    pub fn get_chunk(&self, col: i32, chunk_size: i32) 
    -> PyResult<Vec<String>> {
        check_empty_alignment(self)?;
        check_col_index(self, col as usize)?;
        let mut sequences: Vec<String> = Vec::new();
        let col = col as usize;
        let chunk_size = chunk_size as usize;
        for seq in self.data.iter() {
            let seq_vec: Vec<char> = seq.chars().collect();
            let sequence: String = seq_vec[col..col+chunk_size].to_vec()        
                .into_iter().collect();
            sequences.push(sequence);
        }
        Ok(sequences)
    }

    pub fn get_chunks(&self, cols: Vec<i32>, chunk_size: i32) 
    -> PyResult<Vec<Vec<String>>> {
        check_empty_alignment(self)?;
        if let Some(x) = cols.iter().max() {
            check_col_index(self, *x as usize)?;
        }
        let mut sequences_vec: Vec<Vec<String>> = Vec::new();
        let chunk_size = chunk_size as usize;
        for col in cols.into_iter().map(|x| x as usize) {
            let mut sequences: Vec<String> = Vec::new();
            for seq in self.data.iter() {
                let seq_vec: Vec<char> = seq.chars().collect();
                let sequence: String = seq_vec[col..col+chunk_size].to_vec()
                    .into_iter().collect();
                sequences.push(sequence);
            }
            sequences_vec.push(sequences);
        }
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
    fn remove_cols(&mut self, mut cols: Vec<i32>) -> PyResult<()> {
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
    fn reorder_cols(&mut self, cols: Vec<i32>) -> PyResult<()> {
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


    fn concat(&mut self, others: Vec<&BaseAlignment>) -> PyResult<()> {
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
        Ok(BaseAlignment{ data: self.data.clone() })
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
    m.add_class::<BaseAlignment>()?;

    Ok(())
}
