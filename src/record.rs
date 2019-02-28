use pyo3::prelude::*;
use pyo3::PyObjectProtocol;


#[pyclass(subclass)]
#[derive(Clone)]
/// BaseRecord(id, description, sequence)
/// 
/// BaseRecord represents a single sequence sample.
pub struct BaseRecord {
    #[prop(get,set)]
    pub id: String,

    #[prop(get,set)]
    pub description: String,
    
    #[prop(get,set)]
    pub sequence: String,

}

#[pymethods]
impl BaseRecord {
    #[new]
    /// Creates a new BaseRecord object from id, sequence_description
    /// and sequence.
    fn __new__(obj: &PyRawObject, id: &str, description: &str,
               sequence: &str) -> PyResult<()> {
        obj.init(|_| {
            BaseRecord {
                id: id.to_string(),
                description: description.to_string(),
                sequence: sequence.to_string(),
            }
        })
    }

    #[getter]
    pub fn len(&self) -> PyResult<i32> {
        Ok(self.sequence.len() as i32)
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseRecord {
    fn __repr__(&self) -> PyResult<String> {
        // threshold is 15, 6 ... 6
        let seq: String = match self.len() {
            Ok(x) if x >= 15 => {
                let mut seq = String::new();
                for (i, c) in self.sequence.char_indices() {
                    if i < 6 {
                        seq.push(c);
                    } else if i >= 6 && i < 6 + 3 {
                        seq.push('.');
                    } else if i >= 15 - 6 && i < 15 {
                        seq.push(c);
                    }
                }
                seq
            },
            _ => self.sequence()?.clone()
        };
        Ok(format!(
            "Record(id={id}, sequence=\"{seq}\", length={len})",
            id=self.id, seq=seq, len=self.len()?
        ))
    }

    // default is to output as fasta format
    fn __str__(&self) -> PyResult<String> {
        if self.description.len() > 0 {
            return Ok(format!(">{id} {desc}\n{seq}",
                id=self.id,
                desc=self.description,
                seq=self.sequence))
        }
        return Ok(format!(">{id}\n{seq}",
                id=self.id,
                seq=self.sequence))
    }
}

// Register python functions to PyO3
#[pymodinit]
fn record(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BaseRecord>()?;

    Ok(())
}
