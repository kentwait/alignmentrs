use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};
// use ndarray::{ArrayD, ArrayViewD, ArrayViewMutD};
// use numpy::{IntoPyArray, PyArrayDyn};

use std::str;
use std::io::Write;

#[pyclass]
#[derive(Copy, Clone)]
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
        obj.init(|_| {
            AlignmentMatrix { rows, cols }
        })
    }

}
