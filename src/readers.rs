use pyo3::prelude::*;
use pyo3::exceptions;

use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::Regex;

use crate::record::Record;
use crate::alignment::BaseAlignment;


lazy_static! {
    static ref WS: Regex = Regex::new(r"\s+").unwrap();
}

// FASTA file readers

#[pyfunction]
/// fasta_file_to_records(data_str, /)
/// --
/// 
/// Reads FASTA file and creates a list of Record objects.
fn fasta_file_to_records(path: &str)
-> PyResult<(Vec<Record>, Vec<String>)> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let f = match File::open(path) {
        Err(x) => return Err(exceptions::IOError::py_err(format!("encountered an error while trying to open file {:?}: {:?}", path, x.kind()))),
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
        let sequence = sequence.clone();
        records.push(Record { id, description, sequence });
    }
    Ok((records, comments))
}

#[pyfunction]
/// fasta_file_to_basealignment(data_str, /)
/// --
/// 
/// Reads FASTA file and create a BaseAlignment object.
fn fasta_file_to_basealignment(path: &str)
-> PyResult<(BaseAlignment, Vec<String>)> {  // TODO: Add comments to BaseAlignment structure
    match fasta_file_to_basealignments(path, "") {
        Ok((aln, _, comments)) => Ok((aln, comments)),
        Err(x) => Err(x)
    }
}

#[pyfunction]
/// fasta_file_to_basealignments(data_str, /)
/// --
/// 
/// Reads FASTA file and creates marker and sequence BaseAlignments.
fn fasta_file_to_basealignments(path: &str, marker_kw: &str) -> 
        PyResult<(BaseAlignment, BaseAlignment, Vec<String>)> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let f = match File::open(path) {
        Err(x) => return Err(exceptions::IOError::py_err(
            format!("encountered an error while trying to open file {:?}: {:?}",
                    path, x.kind()))),
        Ok(x) => x
    };
    let f = BufReader::new(f);

    // Declare variables
    let mut s_ids: Vec<String> = Vec::new();
    let mut s_descriptions: Vec<String> = Vec::new();
    let mut s_sequences: Vec<String> = Vec::new();

    let mut m_ids: Vec<String> = Vec::new();
    let mut m_descriptions: Vec<String> = Vec::new();
    let mut m_sequences: Vec<String> = Vec::new();

    let mut comments: Vec<String> = Vec::new();

    // Declare temp variables
    let mut id = String::new();
    let mut description = String::new();
    let mut sequence = String::new();

    // Match regexp
    for line in f.lines() {
        let line = match line {
            Err(x) => return Err(exceptions::IOError::py_err(
                format!("encountered an error while reading file {:?}: {:?}",
                        path, x.kind()))),
            Ok(x) => x.trim().to_string()
        };
        // Handle identifier line \>
        if line.starts_with(">") {
            if sequence.len() > 0 {
                if marker_kw != "" && id.contains(marker_kw) {
                    m_ids.push(id.clone());
                    m_descriptions.push(description.clone());
                    m_sequences.push(sequence.clone());
                } else {
                    s_ids.push(id.clone());
                    s_descriptions.push(description.clone());
                    s_sequences.push(sequence.clone());
                }
                sequence.clear();
            }
            let matches: Vec<&str> = WS.splitn(
                line.trim_start_matches(">"), 2).collect();
            id = matches[0].to_string();
            description = match matches.len() {
                l if l == 2 => matches[1].to_string(),
                _ => String::new(),
            };
        // Handle comment line \;
        } else if line.starts_with(";") {
            id = line.trim_start_matches(";").to_string();
            comments.push(line);
        } else {
            sequence.push_str(&line);
        }
    }
    // Append last sequence line
    if sequence.len() > 0 {
        if marker_kw != "" && id.contains(marker_kw) {
            m_ids.push(id);
            m_descriptions.push(description);
            m_sequences.push(sequence.clone());
        } else {
            s_ids.push(id);
            s_descriptions.push(description);
            s_sequences.push(sequence.clone());
        }
        sequence.clear();
    }
    let sample_aln = BaseAlignment {
        ids: s_ids,
        descriptions: s_descriptions,
        sequences: s_sequences,
    };
    let marker_aln = BaseAlignment {
        ids: m_ids,
        descriptions: m_descriptions,
        sequences: m_sequences,
    };
    Ok((sample_aln, marker_aln, comments))
}

// TODO: Make readers for other file types: PHYLIP, NEXUS


// Register python functions to PyO3
#[pymodinit]
fn readers(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(fasta_file_to_records))?;
    m.add_function(wrap_function!(fasta_file_to_basealignment))?;
    m.add_function(wrap_function!(fasta_file_to_basealignments))?;

    Ok(())
}