use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use ndarray::{ArrayD, ArrayViewD, ArrayViewMutD};
// use numpy::{IntoPyArray, PyArrayDyn};

use std::str;
use std::io::Write;

#[pyclass]
// #[derive(Copy, Clone)]
/// AlignmentMatrix(rows, cols)
/// 
/// Sequence represents a single biological sequence.
pub struct AlignmentMatrix {

    #[prop(get, set)]
    pub rows: usize,

    #[prop(get, set)]
    pub cols: usize,
    
    pub uint_data: Vec<u32>,

}

#[pymethods]
impl AlignmentMatrix {
    #[new]
    /// Creates a new AlignmentMatrix with shape (rows, cols).
    fn __new__(obj: &PyRawObject, rows: usize, cols: usize) -> PyResult<()> {
        let uint_data: Vec<u32> = vec![0; rows * cols];
        obj.init(|_| {
            AlignmentMatrix { rows, cols, uint_data }
        })
    }

}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for AlignmentMatrix {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("AlignmentMatrix(rows={rows}, cols={cols})", rows=self.rows, cols=self.cols))
    }

    // fn __str__(&self) -> PyResult<String> {
    //     Ok(())
    // }
}

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {

    // Add Block class
    m.add_class::<AlignmentMatrix>()?;

    Ok(())
}