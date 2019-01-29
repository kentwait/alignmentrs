use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

use std::str;
use std::io::Write;

#[pyclass]
// #[derive(Copy, Clone)]
/// Sequence(id, description, sequence_str)
/// 
/// Sequence represents a single biological sequence.
pub struct Sequence {

    #[prop(get, set)]
    pub id: String,

    #[prop(get, set)]
    pub description: String,
    
    #[prop(get)]
    pub sequence: String,

}

#[pymethods]
impl Sequence {
    #[new]
    /// Creates a new Sequence object from sequence_id, sequence_description
    /// and sequence_str.
    fn __new__(obj: &PyRawObject, id: &str, description: &str, sequence_str: &str) -> PyResult<()> {
        let sequence_str = String::from(sequence_str);
        obj.init(|_| {
            Sequence {
                id: id.to_string(),
                description: description.to_string(),
                sequence: sequence_str.to_string(),
            }
        })
    }
    fn sequence_to_uint32(&self) -> Vec<i32> {
        let uints: Vec<i32> = self.sequence.chars().map(|x| {x as usize} as i32 ).collect();
        uints
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for Sequence {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("Sequence(id={id}, description={desc}, len={seq_len})", id=self.id, desc=self.description, seq_len=self.sequence.len()))
    }

    fn __str__(&self) -> PyResult<String> {
        if self.description.len() > 0 {
            return Ok(format!(">{id} {desc}\n{seq_len}",
                id=self.id,
                desc=self.description,
                seq_len=self.sequence.len()))
        }
        return Ok(format!(">{id}\n{seq_len}",
                id=self.id,
                seq_len=self.sequence.len()))        
    }
}

// Register python functions to PyO3
#[pymodinit]
fn sequence(_py: Python, m: &PyModule) -> PyResult<()> {

    // Add Block class
    m.add_class::<Sequence>()?;

    Ok(())
}