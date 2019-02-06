use regex::Regex;
use std::io::{BufReader, Read};

lazy_static! {
    static ref FASTA_REGEX: Regex = Regex::new(r"\x3E(?P<id>\S+)\s(?P<description>[^\n]+)\n(?P<sequence>[A-Za-z0-9\n]+)").unwrap();
}

#[pyfunction]
/// fasta_to_lists(data_str)
/// 
/// Reads FASTA file and creates lists of ids, descriptions, and sequences.
fn fasta_to_lists(path: &str) -> PyResult<(Vec<String>, Vec<String>, Vec<String>)> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let mut fasta_str = String::new();
    let f = match File::open(path) {
        Err(x) => return Err(exceptions::IOError::py_err(format!("encountered an error while trying to open file {:?}: {:?}", path, error))),
        _ => ()
    };
    let mut br = BufReader::new(f);
    match br.read_to_string(&mut data) {
        Err(x) => return Err(exceptions::IOError::py_err(format!("encountered an error while reading file {:?} to string: {:?}", path, error))),
        _ => ()
    };

    // Declare variables
    let mut ids: Vec<String> = Vec::new();
    let mut descriptions: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();

    // Match regexp
    for caps in BLOCK_STR_REGEX.captures_iter(fasta_str) {
        ids.push(caps["id"]);
        descriptions.push(caps["description"]);
        sequences.push(caps["sequence"]);
    }
    Ok((ids, descriptions, sequences))
}

// Register python functions to PyO3
#[pymodinit]
fn alignment(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(fasta_to_lists))?;

    Ok(())
}