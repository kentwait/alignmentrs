use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

use std::str;
use std::io::Write;

#[pyclass]
#[derive(Copy, Clone)]
/// Sequence(sequence_id, sequence_description, sequence_str)
/// 
/// Sequence represents a single biological sequence.
pub struct Sequence {

    #[prop(get, set)]
    pub id: i32,

    #[prop(get, set)]
    pub description: String,
    
    #[prop(get)]
    pub uint_sequence: Vec<u32>,

}

#[pymethods]
impl Sequence {
    #[new]
    /// Creates a new Sequence object from sequence_id, sequence_description
    /// and sequence_str.
    fn __new__(obj: &PyRawObject, sequence_id: &str, sequence_description: &str, sequence_str: &str) -> PyResult<()> {
        obj.init(|_| {
            Sequence { sequence_id, sequence_description, sequence_str }
        })
    }

    #[setter(uint_sequence)]
    fn set_uint_sequence(&mut self, uint_sequence: &Vec<u32>) -> PyResult<()> {
        self.uint_sequence = uint_sequence.clone();
        Ok(())
    }
}
