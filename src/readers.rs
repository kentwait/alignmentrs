use pyo3::prelude::*;
use pyo3::exceptions;

use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::Regex;

use crate::record::Record;

lazy_static! {
    static ref WHITESPACE_REGEX: Regex = Regex::new(r"\s+").unwrap();
}

// FASTA file readers

#[pyfunction]
/// fasta_file_to_records(data_str, /)
/// --
/// 
/// Reads FASTA file and creates a list of Record objects.
fn fasta_to_records(path: &str) -> PyResult<(Vec<Record>, Vec<String>)> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let f = match File::open(path) {
        Err(_) => {
            return Err(exceptions::IOError::py_err(format!(
                "encountered an error while trying to open file {:?}", path)))
        },
        Ok(x) => x
    };
    let f = BufReader::new(f);

    // Declare variables
    let mut records: Vec<Record> = Vec::new();
    let mut comments: Vec<String> = Vec::new();

    // Declare temp variables
    let mut id = String::new();
    let mut description = String::new();
    let mut sequence = String::new();

    // Match regexp
    for line in f.lines() {
        let line = match line {
            Err(_) => {
                return Err(exceptions::IOError::py_err(format!(
                    "encountered an error while reading file {:?}", path)))
            },
            Ok(x) => x.trim().to_string()
        };
        if line.starts_with(">") {
            if sequence.len() > 0 {
                let id = id.clone();
                let description = description.clone();
                let sequence: String = sequence.to_string();
                records.push(Record{ id, description, sequence});
            }
            let matches: Vec<&str> = WHITESPACE_REGEX
                .splitn(line.trim_start_matches(">"), 2)
                .collect();
            id = matches[0].to_string();
            description = match matches.len() {
                l if l == 2 => matches[1].to_string(),
                _ => String::new(),
            };
            sequence.clear();
        // Handle comment line \;
        } else if line.starts_with(";") {
            id = line.trim_start_matches(";").to_string();
            comments.push(line);
        } else {
            sequence.push_str(&line);
        }
    }
    if sequence.len() > 0 {
        let id = id.clone();
        let description = description.clone();
        let sequence: String = sequence.to_string();
        records.push(Record{ id, description, sequence});
    }
    Ok((records, comments))
}

// TODO: Make readers for other file types: PHYLIP, NEXUS


// Register python functions to PyO3
#[pymodinit]
fn readers(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(fasta_to_records))?;

    Ok(())
}