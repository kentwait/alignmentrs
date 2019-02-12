use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::Regex;


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
        Ok(format!("Record(id=\"{id}\", len={seq_len}, description=\"{desc}\")", id=self.id, seq_len=self.sequence.len(), desc=desc))
    }

    fn __str__(&self) -> PyResult<String> {
        if self.description.len() > 0 {
            return Ok(format!(">{id} {desc}\n{seq_len}",
                id=self.id,
                desc=self.description,
                seq_len=self.sequence))
        }
        return Ok(format!(">{id}\n{seq_len}",
                id=self.id,
                seq_len=self.sequence))        
    }
}

lazy_static! {
    static ref WS: Regex = Regex::new(r"\s+").unwrap();
}

#[pyfunction]
/// fasta_file_to_records(data_str)
/// 
/// Reads FASTA file and creates a list of Record objects.
fn fasta_file_to_records(path: &str) -> 
        PyResult<Vec<Record>> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let f = match File::open(path) {
        Err(x) => return Err(exceptions::IOError::py_err(format!("encountered an error while trying to open file {:?}: {:?}", path, x.kind()))),
        Ok(x) => x
    };
    let f = BufReader::new(f);

    // Declare variables
    let mut records: Vec<Record> = Vec::new();

    let mut id = String::new();
    let mut description = String::new();
    let mut sequence = String::new();

    // Match regexp
    for line in f.lines() {
        let line = match line {
            Err(x) => return Err(exceptions::IOError::py_err(format!("encountered an error while reading file {:?}: {:?}", path, x.kind()))),
            Ok(x) => x.trim().to_string()
        };
        if line.starts_with(">") {
            if sequence.len() > 0 {
                let id = id.clone();
                let description = description.clone();
                let sequence = sequence.clone();
                records.push(Record { id, description, sequence });
            }
            let matches: Vec<&str> = WS.splitn(line.trim_start_matches(">"), 2).collect();
            id = matches[0].to_string();
            description = match matches.len() {
                l if l == 2 => matches[1].to_string(),
                _ => String::new(),
            };
            sequence.clear();
        } else {
            sequence.push_str(&line);
        }
    }
    if sequence.len() > 0 {
        let id = id.clone();
        let description = description.clone();
        let sequence = sequence.clone();
        records.push(Record { id, description, sequence });
    }
    Ok(records)
}

// Register python functions to PyO3
#[pymodinit]
fn record(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Record>()?;
    m.add_function(wrap_function!(fasta_file_to_records))?;

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::{PyObjectProtocol, Record};

    #[test]
    fn test_repr() {
        let expected = "Record(id=\"test\", len=9, description=\"Test description\")";
        let record = Record{
            id: "test".to_string(),
            description: "Test description".to_string(),
            sequence: "ATGCGATTA".to_string()
        };
        let actual = record.__repr__().expect("__repr__ method failed");
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_str() {
        let expected = ">test Test description\nATGCGATTA";
        let record = Record{
            id: "test".to_string(),
            description: "Test description".to_string(),
            sequence: "ATGCGATTA".to_string()
        };
        let actual = record.__str__().expect("__str__ method failed");
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_len() {
        let expected = 9;
        let record = Record{
            id: "test".to_string(),
            description: "Test description".to_string(),
            sequence: "ATGCGATTA".to_string()
        };
        let actual = record.len().expect("len method failed");
        assert_eq!(expected, actual);
    }
}