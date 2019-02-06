use pyo3::prelude::*;
use pyo3::{PyObjectProtocol};


#[pyclass]
// #[derive(Copy, Clone)]
/// Record(id, description, sequence_str)
/// 
/// Record represents a single biological sequence.
pub struct Record {

    #[prop(get, set)]
    pub id: String,

    #[prop(get, set)]
    pub description: String,
    
    #[prop(get)]
    pub sequence: String,

}

#[pymethods]
impl Record {
    #[new]
    /// Creates a new Record object from sequence_id, sequence_description
    /// and sequence_str.
    fn __new__(obj: &PyRawObject, id: &str, description: &str, sequence_str: &str) -> PyResult<()> {
        let sequence_str = String::from(sequence_str);
        obj.init(|_| {
            Record {
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
impl PyObjectProtocol for Record {
    fn __repr__(&self) -> PyResult<String> {
        let desc:String = match self.description.chars().count() {
            x if x > 20 => {
                let mut desc:String = self.description.char_indices()
                                        .filter(|(i, _)| *i < 17)
                                        .map(|(_, c)| c)
                                        .collect();
                desc.push_str("...");
                desc
            },
            _ => self.description.clone()
        };
        Ok(format!("Record(id=\"{id}\", len={seq_len}, description=\"{desc}\")", id=self.id, seq_len=self.sequence.len(), desc=desc))
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
fn record(_py: Python, m: &PyModule) -> PyResult<()> {

    // Add Block class
    m.add_class::<Record>()?;

    Ok(())
}