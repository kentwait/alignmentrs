use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};


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
    
    pub sequence: Vec<String>,

    #[prop(get)]
    pub chunk_size: i32,

}

#[pymethods]
impl BaseRecord {
    #[new]
    /// Creates a new BaseRecord object from sequence_id, sequence_description
    /// and sequence_str.
    fn __new__(obj: &PyRawObject, id: &str, description: &str, sequence_str: &str, chunk_size: i32) -> PyResult<()> {
        if sequence_str.chars().count() % (chunk_size as usize) != 0 {
            return Err(exceptions::IndexError::py_err(
                "invalid chunk_size: sequence cannot be cleanly divided into substrings"))
        }
        obj.init(|_| {
            let chars: Vec<char> = sequence_str.chars().collect();
            let sequence: Vec<String> = chars.chunks(chunk_size as usize)
                .map(|chunk| chunk.iter().collect::<String>())
                .collect();
            BaseRecord {
                id: id.to_string(),
                description: description.to_string(),
                sequence: sequence,
                chunk_size: chunk_size,
            }
        })
    }

    #[getter]
    pub fn len(&self) -> PyResult<i32> {
        Ok((self.sequence.len() as i32)  * self.chunk_size)
    }

    #[getter]
    pub fn chunked_len(&self) -> PyResult<i32> {
        return Ok(self.sequence.len() as i32)
    }

    #[setter(chunk_size)]
    pub fn set_chunk_size(&mut self, chunk_size: i32) -> PyResult<()> {
        if self.len()? % chunk_size != 0 {
            return Err(exceptions::IndexError::py_err(
                "invalid chunk_size: sequence cannot be cleanly divided into substrings"))
        }
        let chars: Vec<char> = self.sequence.join("").chars().collect();
        let sequence: Vec<String> = chars.chunks(chunk_size as usize)
            .map(|chunk| chunk.iter().collect::<String>())
            .collect();
        self.sequence = sequence;
        Ok(())
    }

    #[getter(sequence)]
    pub fn get_sequence(&self) -> PyResult<String> {
        Ok(self.sequence.join(""))
    }

    #[setter(sequence)]
    pub fn set_sequence(&mut self, sequence_str: &str) -> PyResult<()> {
        let chars: Vec<char> = sequence_str.chars().collect();
        let sequence: Vec<String> = chars.chunks(self.chunk_size as usize)
            .map(|chunk| chunk.iter().collect::<String>())
            .collect();
        self.sequence = sequence;
        Ok(())
    }

    #[getter]
    pub fn chunked_sequence(&self) -> PyResult<Vec<String>> {
        Ok(self.sequence.clone())
    }

    // TODO: Get/set chunk
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseRecord {
    fn __repr__(&self) -> PyResult<String> {
        // threshold is 15, 6 ... 6
        let seq: String = match self.len() {
            Ok(x) if x >= 15 => {
                let sequence = self.get_sequence()?;
                let mut seq = String::new();
                for (i, c) in sequence.char_indices() {
                    if i < 30 {
                        seq.push(c);
                    } else if i >= 6 && i < 6 + 3 {
                        seq.push('.');
                    } else if i >= 15 - 6 && i < 15 {
                        seq.push(c);
                    }
                }
                seq
            },
            _ => self.get_sequence()?.clone()
        };
        Ok(format!(
            "BaseRecord(id={id}, sequence=\"{seq}\", length={len}, chunk_size={chunk_size}",
            id=self.id, seq=seq, len=self.len()?, chunk_size=self.chunk_size
        ))
    }

    // default is to output as fasta format
    fn __str__(&self) -> PyResult<String> {
        if self.description.len() > 0 {
            return Ok(format!(">{id} {desc}\n{seq_len}",
                id=self.id,
                desc=self.description,
                seq_len=self.len()?))
        }
        return Ok(format!(">{id}\n{seq_len}",
                id=self.id,
                seq_len=self.len()?))
    }
}

// Register python functions to PyO3
#[pymodinit]
fn record(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BaseRecord>()?;

    Ok(())
}
