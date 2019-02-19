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
    
    #[prop(get,set)]
    pub sequence: String,

    #[prop(get,set)]
    pub chunk_size: i32,

}

#[pymethods]
impl BaseRecord {
    #[new]
    /// Creates a new BaseRecord object from sequence_id, sequence_description
    /// and sequence_str.
    fn __new__(obj: &PyRawObject, id: &str, description: &str, sequence_str: &str, chunk_size: i32) -> PyResult<()> {
        let sequence_str = String::from(sequence_str);
        obj.init(|_| {
            BaseRecord {
                id: id.to_string(),
                description: description.to_string(),
                sequence: sequence_str.to_string(),
                chunk_size: chunk_size,
            }
        })
    }

    #[getter]
    fn len(&self) -> PyResult<i32> {
        Ok(self.sequence.chars().count() as i32)
    }

    #[getter]
    fn chunklen(&self) -> PyResult<i32> {
        let nchunks = self.sequence.chars().count() % self.chunk_size as usize;
        if nchunks != 0 {
            return Err(exceptions::ValueError::py_err(
                "sequence cannot be cleanly divided into chunks."))
        }
        return Ok((self.sequence.chars().count() / self.chunk_size as usize) as i32)
    }

    fn chunked(&self) -> PyResult<Vec<String>> {
        if self.sequence.len() == 0 {
            return Ok(Vec::new())
        }
        let nchunks = self.sequence.chars().count() % self.chunk_size as usize;
        if nchunks != 0 {
            return Err(exceptions::ValueError::py_err(
                "sequence cannot be cleanly divided into chunks."))
        }
        let mut chunks: Vec<String> = Vec::with_capacity(nchunks);
        // temp variables and counter
        let mut chunk: String = String::new();
        let mut i = 0;
        // iterate over characters
        for c in self.sequence.chars() {
            if i < self.chunk_size {
                chunk.push(c);
            }
            i += 1;
            // Add chunk to list and reset temp variables and counter
            if i == self.chunk_size {
                chunks.push(chunk);
                chunk.clear();
                i = 0;
            }
        }
        Ok(chunks)
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for BaseRecord {
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
        Ok(format!("BaseRecord(id=\"{id}\", sequence=\"{seq}\", len={seq_len}, \
                   chunk_size={chunk_size}, description=\"{desc}\")", 
                   id=self.id, seq=seq, seq_len=self.sequence.len(),
                   chunk_size=self.chunk_size, desc=desc))
    }

    // default is to output as fasta format
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


// Register python functions to PyO3
#[pymodinit]
fn BaseRecord(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BaseRecord>()?;

    Ok(())
}
