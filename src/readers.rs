use pyo3::prelude::*;
use pyo3::exceptions;

use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;
use regex::Regex;

use crate::alignment::SeqMatrix;

lazy_static! {
    static ref WHITESPACE_REGEX: Regex = Regex::new(r"\s+").unwrap();
}

// FASTA file readers

pub fn fasta_to_hashmap(path: &str) -> Result<HashMap<String, Vec<String>>, String> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let f = match File::open(path) {
        Err(_) => {
            return Err(format!(
                "encountered an error while trying to open file {:?}", path))
        },
        Ok(x) => x
    };
    let f = BufReader::new(f);

    // Declare variables
    let mut ids: Vec<String> = Vec::new();
    let mut descriptions: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut comments: Vec<String> = Vec::new();

    // Declare temp variables
    let mut id = String::new();
    let mut description = String::new();
    let mut sequence = String::new();

    // Match regexp
    for line in f.lines() {
        let line = match line {
            Err(_) => {
                return Err(format!(
                    "encountered an error while reading file {:?}", path))
            },
            Ok(x) => x.trim().to_string()
        };
        if line.starts_with(">") {
            if sequence.len() > 0 {
                ids.push(id.clone());
                descriptions.push(description.clone());
                sequences.push(sequence.clone());
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
        ids.push(id.clone());
        descriptions.push(description.clone());
        sequences.push(sequence.clone());
    }

    let map: HashMap<String, Vec<String>> = [
        ("ids".to_string(), ids),
        ("descriptions".to_string(), descriptions),
        ("sequences".to_string(), sequences),
        ("comments".to_string(), comments),
    ].iter().cloned().collect();

    Ok(map)
}

#[pyfunction]
/// fasta_to_dict(data_str, /)
/// --
/// 
/// Reads FASTA file and creates a list of Record objects.
fn fasta_to_dict(path: &str) -> PyResult<(SeqMatrix, HashMap<String, Vec<String>>)> {
    match fasta_to_hashmap(path) {
        Ok(mut d) => {
            let data = d.remove("sequences").unwrap();
            let rows = data.len();
            let cols = if rows > 0 {data[0].len()} else {0};
            Ok((SeqMatrix{ data, rows, cols }, d))
        },
        Err(x) => Err(exceptions::IOError::py_err(x))
    }
}

// TODO: Make readers for other file types: PHYLIP, NEXUS


// Register python functions to PyO3
#[pymodinit]
fn readers(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(fasta_to_dict))?;

    Ok(())
}