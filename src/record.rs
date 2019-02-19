use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};


#[pyclass(subclass)]
#[derive(Clone)]
/// Record(id, description, sequence)
/// 
/// Record represents a single sequence sample.
pub struct Record {
    #[prop(get,set)]
    pub id: String,

    #[prop(get,set)]
    pub description: String,
    
    #[prop(get,set)]
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

    #[getter]
    fn len(&self) -> PyResult<i32> {
        Ok(self.sequence.chars().count() as i32)
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
        let seq: String = match self.sequence.chars().count() {
            x if x >= 15 => {
                let mut seq = String::new();
                for (i, c) in self.sequence.char_indices() {
                    if i < 6 || i >= x-6 {
                        seq.push(c);
                    } else if i >= 6 && i < 9 {
                        seq.push('.');
                    }
                }
                seq
            },
            _ => self.sequence.clone()
        };
        Ok(format!("Record(id=\"{id}\", sequence=\"{seq}\", len={seq_len}, \
                   description=\"{desc}\")", 
                   id=self.id, seq=seq,
                   seq_len=self.sequence.len(), desc=desc))
    }

    fn __str__(&self) -> PyResult<String> {
        if self.description.len() > 0 {
            return Ok(format!(">{id} {desc}\n{seq_len}",
                id=self.id,
                desc=self.description,
                seq_len=self.sequence.chars().count()))
        }
        return Ok(format!(">{id}\n{seq_len}",
                id=self.id,
                seq_len=self.sequence.chars().count()))        
    }
}

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
    fn len(&self) -> PyResult<i32> {
        Ok((self.sequence.len() as i32)  * self.chunk_size)
    }

    #[getter]
    fn chunked_len(&self) -> PyResult<i32> {
        return Ok(self.sequence.len() as i32)
    }

    #[setter(chunk_size)]
    fn set_chunk_size(&mut self, chunk_size: i32) -> PyResult<()> {
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
    fn get_sequence(&self) -> PyResult<String> {
        Ok(self.sequence.join(""))
    }

    #[setter(sequence)]
    fn set_sequence(&mut self, sequence_str: &str) -> PyResult<()> {
        let chars: Vec<char> = sequence_str.chars().collect();
        let sequence: Vec<String> = chars.chunks(self.chunk_size as usize)
            .map(|chunk| chunk.iter().collect::<String>())
            .collect();
        self.sequence = sequence;
        Ok(())
    }

    #[getter]
    fn chunked_sequence(&self) -> PyResult<Vec<String>> {
        Ok(self.sequence.clone())
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseRecord {
    fn __repr__(&self) -> PyResult<String> {
        // threshold cols is 80
        let desc:String = match self.description.chars().count() {
            // 80 - 13
            x if x > 67 => {
                let mut desc:String = self.description.char_indices()
                                        .filter(|(i, _)| *i < 67 - 3)
                                        .map(|(_, c)| c)
                                        .collect();
                desc.push_str("...");
                desc
            },
            _ => self.description.clone()
        };
        let seq: String = match self.len() {
            // 80 - 10
            Ok(x) if x >= 70 => {
                let sequence = self.get_sequence()?;
                let mut seq = String::new();
                for (i, c) in sequence.char_indices() {
                    if i < 30 {
                        seq.push(c);
                    } else if i >= 30 && i < 33 {
                        seq.push('.');
                    } else if i >= 70 - 30 && i < 70 {
                        seq.push(c);
                    }
                }
                seq
            },
            _ => self.get_sequence()?.clone()
        };
        Ok(format!(
            "[BaseRecord]\n\
            id = {id}\n\
            description = \"{desc}\"\n\
            sequence = \"{seq}\"\n\
            length = {len}\n\
            chunk_size = {chunk_size}\n\
            chunked_len = {chunked_len}",
            id=self.id, desc=desc, seq=seq, len=self.len()?, chunk_size=self.chunk_size,
            chunked_len=self.chunked_len()?
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
    m.add_class::<Record>()?;
    m.add_class::<BaseRecord>()?;

    Ok(())
}
