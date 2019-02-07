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
/// LSpace(init_state, start, stop)
/// 
/// LSpace represents a dimension in space.
pub struct LSpace {

    // type_map: HashMap<i32, String>,
    blocks: Vec<[i32; 3]>,

}

#[pymethods]
impl LSpace {
    #[new]
    /// Creates a new LSpace object from an init_state, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, init_state: i32, start: i32, stop: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            LSpace { 
                blocks: vec![[init_state, start, stop]],
            }
        })
    }

    /// Removes points in linear space given based on a list of coordinates.
    fn remove(&mut self, coords: Vec<i32>) -> PyResult<()> {
        let mut coords = coords;
        coords.sort_unstable();
        coords.reverse();
        let mut offset = 0;
        for pos in coords.iter() {
            let pos = *pos;
            let length = self.blocks.len() - offset;
            for j in (0..length).rev() {
                if pos > self.blocks[j][1] && pos < self.blocks[j][2] - 1 {
                    // Remove block currently at j and get values
                    let [id, start, stop] = self.blocks[j];
                    // Split this block at pos
                    // 10,11,12,13,14,15 : remove at index 2 (3rd pos)
                    // [_,10,16]
                    // 10,11,   13,14,15
                    // [_,10,12], [13:16]
                    // Insert two new blocks at j and j+1
                    self.blocks[j] = [id, start, pos];
                    self.blocks.insert(j+1, [id, pos+1, stop]);
                    offset += 1;
                } else if pos == self.blocks[j][1] {
                    if self.blocks[j][1] == self.blocks[j][2] - 1 {
                        let _ = self.blocks.remove(j);
                        // no offset increment
                    } else {
                        let [id, _, stop] = self.blocks[j];
                        self.blocks[j] = [id, pos+1, stop];
                        // no offset increment
                    }
                } else if pos == self.blocks[j][2] - 1 {
                    let [id, start, _] = self.blocks[j];
                    self.blocks[j] = [id, start, pos];
                    // no offset increment
                }
            }
        }
        Ok(())
    }

    /// Retains points in linear space specified by a
    /// list of coordinates to keep.
    fn retain(&mut self, coords: Vec<i32>) -> PyResult<()> {
        if let Some([_, _, stop]) = self.blocks.last() {
            let inverse_ilist: Vec<i32> = (0..*stop)
                                            .filter(|x| !coords.contains(x))
                                            .collect();
            return self.remove(inverse_ilist)
        };
        Err(exceptions::ValueError::py_err("cannot perform retain on \
                                            dimension: block list is empty"))
    }

    /// Inserts into the linear space at the given position.
    fn insert(&mut self, pos: i32, state: i32, length: i32) -> PyResult<()> {
        for i in 0..self.blocks.len() {
            let [id, start, stop] = self.blocks[i];
            if start >= pos && stop < pos {
                if id == state {
                    // Extend current block
                    self.blocks[i] = [id, start, stop + length];
                } else {
                    // Create new state
                    self.blocks.insert(i + 1, [id, stop, stop + length]);
                }
                // Adjust
                for j in i..self.blocks.len() {
                    let [id, start, stop] = self.blocks[j];
                    self.blocks[j] = [id, start+length, stop+length];
                }
                break;
            }
        }
        Ok(())
    }

    /// Appends to the end of the linear space.
    fn append(&mut self, state: i32, length: i32) -> PyResult<()> {
        if let Some([id, start, stop]) = self.blocks.last() {
            let (id, start, stop) = (*id, *start, *stop);
            let i = self.blocks.len();
            if id == state {
                // Extend current block
                self.blocks[i] = [id, start, stop + length];
            } else {
                // Create new state
                self.blocks.push([id, stop, stop + length]);
            }
        } else {
            self.blocks.push([state, 0, length]);
        }
        Ok(())
    }

    // start, stop, full_len

    /// Returns the value at the start of the linear space.
    fn start(&self) -> PyResult<i32> {
        match self.blocks.first() {
            Some([_, x, _]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err("cannot get minimum position: block list is empty"))
        }
    }

    /// Returns the end value of the linear space.
    /// This value is not part of the space.
    fn stop(&self) -> PyResult<i32> {
        match self.blocks.last() {
            Some([_, _, x]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err("cannot get maximum position: block list is empty"))
        }
    }

    /// Returns the total length of the linear space.
    fn len_all(&self) -> PyResult<i32> {
        let mut length = 0;
        for [_, start, stop] in self.blocks.iter() {
            let (start, stop) = (*start, *stop);
            length += stop - start;
        }
        Ok(length)
    }

    /// Returns the total length of the linear space where the
    /// state is equal to 1.
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

    /// Returns the total length of the linear space where the
    /// state is equal to 0.
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

    #[staticmethod]
    /// Returns a linear space created using the given list of blocks.
    fn from_blocks(blocks: Vec<&Block>) -> PyResult<LSpace> {
        let mut data: Vec<[i32; 3]> = Vec::new();
        for Block{ id, start, stop} in blocks.iter() {
            let id: i32 = match id.parse() {
                Ok(x) => x,
                Err(_) => return Err(exceptions::TypeError::py_err("Cannot parse block id, maybe not i32?"))
            };
            data.push([id, *start, *stop]);
        }
        Ok(LSpace { blocks: data })
    }

    /// Returns the linear space as a list of blocks.
    fn to_blocks(&self) -> PyResult<Vec<Block>> {
        Ok(
            self.blocks.iter()
                .map(|[id, start, stop]| {
                    let id = match *id {
                        1 => "s",
                        _ => "g"
                    }.to_string();
                    Block { id, start: *start, stop: *stop }
                }).collect()
        )
    }

    // TODO: from_array

    /// Returns the linear space as a list of integer coordinates.
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

    /// Returns a deep copy of the current linear space.
    fn copy(&self) -> PyResult<LSpace> {
        let blocks = self.blocks.clone();
        Ok(LSpace{ blocks })
    }
}

#[pyproto]
impl PyObjectProtocol for LSpace {
    fn __repr__(&self) -> PyResult<String> {
        let start = match self.start() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let stop = match self.stop() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let length = match self.len_all() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        Ok(format!("LSpace(start={}, stop={}, length={})", start, stop, length))
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

impl LSpace {
    fn _to_array(&self) -> Vec<i32> {
        let mut coords: Vec<i32> = Vec::new();
        for [id, start, stop] in self.blocks.iter() {
            let (id, start, stop) = (*id, *start, *stop);
            if id == 1 {
                let mut array: Vec<i32> = (start..stop).map(|x| x).collect();
                coords.append(&mut array);
            } else if id == 0 {
                let mut array: Vec<i32> = (start..stop).map(|_| -1).collect();
                coords.append(&mut array);
            }
        }
        coords
    }
}


#[pymodinit]
fn position(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Block>()?;
    m.add_class::<LSpace>()?;

    Ok(())
}