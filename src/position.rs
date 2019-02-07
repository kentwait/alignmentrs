use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

#[pyclass(subclass)]
#[derive(Clone)]
/// Block(id, start, stop)
/// 
/// Block represents a linear interval.
pub struct Block {

    #[prop(get,set)]
    pub id: String,

    #[prop(get,set)]
    pub start: i32,
    
    #[prop(get,set)]
    pub stop: i32,

}

#[pymethods]
impl Block {
    #[new]
    /// Creates a new Block object from an id, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, id: &str, start: i32, stop: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            let id = id.to_string();
            Block { id, start, stop }
        })
    }

    /// Returns True if a given position is inside the block.
    /// Otherwise, returns False.
    fn in_block(&self, i: i32) -> PyResult<bool> {
        if i >= self.start && i < self.stop {
            return Ok(true)
        }
        Ok(false)
    }

    /// Converts the block into a list of positions.
    fn to_array(&self) -> PyResult<Vec<i32>> {
        Ok(self._to_array())
    }

    // Formatting methods

    /// Converts block into a compressed string representation.
    fn to_compressed_str(&self) -> PyResult<String> {
        self.__str__()
    }

    /// Converts block into an extended (human-readable) string
    /// representation.
    fn to_extended_str(&self) -> PyResult<String> {
        Ok(format!("{}={}:{}", self.id, self.start, self.stop))
    }

    /// Converts block into comma-separated list of positions.
    fn to_array_str(&self) -> PyResult<String> {
        let v: Vec<String> = self._to_array().iter()
                                .map(|x| format!("{}", x))
                                .collect();
        Ok(v.join(","))
    }
}

impl Block {
    fn _to_array(&self) -> Vec<i32> {
        (self.start..self.stop).collect::<Vec<i32>>()
    }
}

#[pyproto]
impl PyObjectProtocol for Block {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("Block(id=\"{}\", start={}, stop={}",
                   self.id, self.start, self.stop))
    }
    
    fn __str__(&self) -> PyResult<String> {
        Ok(format!("{}{}{}", self.start, self.id, self.stop))
    }
}

#[pyclass(subclass)]
#[derive(Clone)]
pub struct Dimension {

    // type_map: HashMap<i32, String>,
    blocks: Vec<[i32; 3]>,

}

#[pymethods]
impl Dimension {
    #[new]
    /// Creates a new Dimension object from an id, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, init_type: i32, start: i32, stop: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            Dimension { 
                blocks: vec![[init_type, start, stop]],
            }
        })
    }

    fn remove(&mut self, ilist: Vec<i32>) -> PyResult<()> {
        let mut ilist = ilist;
        ilist.sort_unstable();
        ilist.reverse();
        for (i, pos) in ilist.iter().enumerate() {
            let pos = *pos;
            for j in (0..(self.blocks.len() - i)).rev() {
                if pos >= self.blocks[j][1] && pos < self.blocks[j][2] {
                    // Remove block currently at j and get values
                    let [t, start, stop] = self.blocks.remove(j);
                    // Split this block at pos
                    // 10,11,12,13,14,15 : remove at index 2 (3rd pos)
                    // [_,10,16]
                    // 10,11,   13,14,15
                    // [_,10,12], [13:16]
                    // Insert two new blocks at j and j+1
                    self.blocks.insert(j, [t, start, pos]);
                    self.blocks.insert(j+1, [t, pos+1, stop]);
                }
            }
        }
        Ok(())
    }

    fn retain(&mut self, ilist: Vec<i32>) -> PyResult<()> {
        if let Some([_, _, last]) = self.blocks.last() {
            let inverse_ilist: Vec<i32> = (0..*last)
                                            .filter(|x| !ilist.contains(x))
                                            .collect();
            return self.remove(inverse_ilist)
        };
        Err(exceptions::ValueError::py_err("cannot perform retain on \
                                            dimension: block list is empty"))
    }

    // min, max, full_len
    fn min(&self) -> PyResult<i32> {
        match self.blocks.first() {
            Some([_, _, x]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err("cannot get minimum position: block list is empty"))
        }
    }

    fn max(&self) -> PyResult<i32> {
        match self.blocks.last() {
            Some([_, _, x]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err("cannot get maximum position: block list is empty"))
        }
    }

    fn len_all(&self) -> PyResult<i32> {
        let mut length = 0;
        for [_, start, stop] in self.blocks.iter() {
            let (start, stop) = (*start, *stop);
            length += stop - start;
        }
        Ok(length)
    }

    fn len_seq(&self) -> PyResult<i32> {
        let mut length = 0;
        for [id, start, stop] in self.blocks.iter() {
            let (id, start, stop) = (*id, *start, *stop);
            if id == 1 {
                length += stop - start;
            }
        }
        Ok(length)
    }

    fn len_gap(&self) -> PyResult<i32> {
        let mut length = 0;
        for [id, start, stop] in self.blocks.iter() {
            let (id, start, stop) = (*id, *start, *stop);
            if id == 0 {
                length += stop - start;
            }
        }
        Ok(length)
    }
    

    // Format conversion

    fn to_blocks(&self) -> PyResult<Vec<Block>> {
        Ok(
            self.blocks.iter()
                .map(|[id, start, stop]| {
                    Block { id: format!("{}", id), start: *start, stop: *stop }
                }).collect()
        )
    }

    fn to_array(&self) -> PyResult<Vec<i32>> {
        Ok(self._to_array())
    }

    // Formatting methods

    /// Converts block into a compressed string representation.
    fn to_compressed_str(&self) -> PyResult<String> {
        self.__str__()
    }

    /// Converts block into an extended (human-readable) string
    /// representation.
    fn to_extended_str(&self) -> PyResult<String> {
        let mut strings: Vec<String> = Vec::new();
        if let Ok(blocks) = self.to_blocks() {
            for block in blocks {
                if let Ok(s) = block.to_extended_str() {
                    strings.push(s);
                } else {
                    return Err(exceptions::ValueError::py_err("cannot get string representation of block"))
                }
            }
        } else {
            return Err(exceptions::ValueError::py_err("cannot generate blocks"))
        }
        Ok(strings.join(","))
    }

    /// Converts block into comma-separated list of positions.
    fn to_array_str(&self) -> PyResult<String> {
        let mut strings: Vec<String> = Vec::new();
        if let Ok(blocks) = self.to_blocks() {
            for block in blocks {
                if let Ok(s) = block.to_array_str() {
                    strings.push(s);
                } else {
                    return Err(exceptions::ValueError::py_err("cannot get string representation of block"))
                }
            }
        } else {
            return Err(exceptions::ValueError::py_err("cannot generate blocks"))
        }
        Ok(strings.join(","))
    }
}

#[pyproto]
impl PyObjectProtocol for Dimension {
    fn __repr__(&self) -> PyResult<String> {
        let min = match self.min() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let max = match self.max() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let length = match self.len_all() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        Ok(format!("Dimension(min={}, max={}, length={})", min, max, length))
    }
    
    fn __str__(&self) -> PyResult<String> {
        let mut strings: Vec<String> = Vec::new();
        if let Ok(blocks) = self.to_blocks() {
            for block in blocks {
                if let Ok(s) = block.__str__() {
                    strings.push(s);
                } else {
                    return Err(exceptions::ValueError::py_err("cannot get string representation of block"))
                }
            }
        } else {
            return Err(exceptions::ValueError::py_err("cannot generate blocks"))
        }
        Ok(strings.join(","))
    }
}

impl Dimension {
    fn _to_array(&self) -> Vec<i32> {
        let mut ilist: Vec<i32> = Vec::new();
        for [id, start, stop] in self.blocks.iter() {
            let (id, start, stop) = (*id, *start, *stop);
            if id == 1 {
                let mut array: Vec<i32> = (start..stop).map(|x| x).collect();
                ilist.append(&mut array);
            } else if id == 0 {
                let mut array: Vec<i32> = (start..stop).map(|_| -1).collect();
                ilist.append(&mut array);
            }
        }
        ilist
    }
}


#[pymodinit]
fn position(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Block>()?;
    m.add_class::<Dimension>()?;

    Ok(())
}