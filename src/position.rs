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
/// BlockSpace(init_state, start, stop)
/// 
/// BlockSpace represents a discrete linear space stored
/// in intervals called Blocks.
pub struct BlockSpace {

    // type_map: HashMap<i32, String>,
    blocks: Vec<[i32; 3]>,

}

#[pymethods]
impl BlockSpace {
    #[new]
    /// Creates a new BlockSpace object from an init_state, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, init_state: i32, start: i32, stop: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            BlockSpace { 
                blocks: vec![[init_state, start, stop]],
            }
        })
    }

    // These methods relies on reading into the block
    // But another way to access it is the relative position of the block

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
    fn from_arrays(data: Vec<i32>, ids: Vec<String>) -> PyResult<BlockSpace> {
        match arrays_to_blocks(data, ids) {
            Ok(x) => {
                let mut blocks: Vec<[i32; 3]> = Vec::new();
                for Block{ id, start, stop } in x.iter() {
                    let start = *start;
                    let stop = *stop;
                    if id == "s" {
                        blocks.push([1, start, stop]);
                    } else if id == "g" {
                        blocks.push([0, start, stop]);
                    } else {
                        return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", id)))
                    }
                }
                Ok(BlockSpace{ blocks })
            },
            Err(x) => Err(x),
        }
    }

    #[staticmethod]
    /// Returns a linear space created using the given list of blocks.
    fn from_blocks(blocks: Vec<&Block>) -> PyResult<BlockSpace> {
        let mut data: Vec<[i32; 3]> = Vec::new();
        for Block{ id, start, stop} in blocks.iter() {
            let id: i32 = match id.parse() {
                Ok(x) => x,
                Err(_) => return Err(exceptions::TypeError::py_err("Cannot parse block id, maybe not i32?"))
            };
            data.push([id, *start, *stop]);
        }
        Ok(BlockSpace { blocks: data })
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

    /// Returns the linear space as a list of integer coordinates.
    fn to_arrays(&self) -> PyResult<(Vec<i32>, Vec<String>)> {
        let mut data: Vec<i32> = Vec::new();
        let mut ids: Vec<String> = Vec::new();
        for [id, start, stop] in self.blocks.iter() {
            let (id, start, stop) = (*id, *start, *stop);
            if id == 1 {
                let mut data_array: Vec<i32> = (start..stop).map(|x| x).collect();
                let mut id_array: Vec<String> = vec!["s".to_string(); data_array.len()];
                data.append(&mut data_array);
                ids.append(&mut id_array);
            } else if id == 0 {
                let mut data_array: Vec<i32> = (start..stop).map(|_| -1).collect();
                let mut id_array: Vec<String> = vec!["g".to_string(); data_array.len()];
                data.append(&mut data_array);
                ids.append(&mut id_array);
            } else {
                return Err(exceptions::ValueError::py_err(format!("unexpected id value: {}", id)))
            }

        }
        Ok((data, ids))
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
    fn copy(&self) -> PyResult<BlockSpace> {
        let blocks = self.blocks.clone();
        Ok(BlockSpace{ blocks })
    }
}

#[pyproto]
impl PyObjectProtocol for BlockSpace {
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
        Ok(format!("BlockSpace(start={}, stop={}, length={})", start, stop, length))
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


#[pyclass(subclass)]
#[derive(Clone)]
/// CoordSpace(init_state, start, stop)
/// 
/// CoordSpace represents a discrete linear space stored as a 
/// list of integer coordinates.
pub struct CoordSpace {

    coords: Vec<i32>

}

#[pymethods]
impl CoordSpace {
    #[new]
    /// Creates a new CoordSpace object from an init_state, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, start: i32, stop: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            CoordSpace { 
                coords: (start..stop).collect(),
            }
        })
    }

    /// Removes points in linear space given based on a list of relative
    /// coordinates.
    fn remove(&mut self, coords: Vec<i32>) -> PyResult<()> {
        self.coords = self.coords.iter().enumerate().filter(|(i, _)| coords.contains(&(*i as i32))).map(|(_, x)| *x ).collect();
        Ok(())
    }


    /// Retains points in linear space specified by a
    /// list of coordinates to keep.
    fn retain(&mut self, coords: Vec<i32>) -> PyResult<()> {
        self.coords = self.coords.iter().enumerate().filter(|(i, _)| !coords.contains(&(*i as i32))).map(|(_, x)| *x ).collect();
        Ok(())
    }

    /// Inserts into the linear space at the given position.
    fn insert(&mut self, pos: i32, start: i32, length: i32) -> PyResult<()> {
        // Insert to start of list if pos is 0,
        // Append to end of list if pos is the length,
        // Otherwise split at the position and concat
        if pos == 0 {
            let mut new_coords: Vec<i32> = (start..start+length).collect();
            new_coords.append(&mut self.coords);
            self.coords = new_coords;
            return Ok(())
        } else if pos == self.coords.len() as i32 {
            return self.append(start, length)
        }
        // Split list into two and combine
        let mut new_coords: Vec<i32> = (start..start+length).collect();
        let (left, right) = self.coords.split_at(pos as usize);
        let mut left: Vec<i32> = left.iter().map(|x| *x ).collect();
        let mut right: Vec<i32> = right.iter().map(|x| *x ).collect();
        left.append(&mut new_coords);
        left.append(&mut right);
        self.coords = left;
        Ok(())
    }

    /// Appends to the end of the linear space.
    fn append(&mut self, start: i32, length: i32) -> PyResult<()> {
        let mut new_coords: Vec<i32> = (start..start+length).collect();
        self.coords.append(&mut new_coords);
        Ok(())
    }

    // start, stop, full_len

    /// Returns the value at the start of the linear space.
    fn start(&self) -> PyResult<i32> {
        Ok(self.coords[0])
    }

    /// Returns the end value of the linear space.
    /// This value is not part of the space.
    fn stop(&self) -> PyResult<i32> {
        Ok(self.coords[self.coords.len()])
    }

    /// Returns the total length of the linear space.
    fn len_all(&self) -> PyResult<i32> {
        Ok(self.coords.len() as i32)
    }

    /// Returns the total length of the linear space where the
    /// state is equal to 1.
    fn len_seq(&self) -> PyResult<i32> {
        let length = self.coords.iter().filter(|x| **x > 0).collect::<Vec<&i32>>().len();
        Ok(length as i32)
    }

    /// Returns the total length of the linear space where the
    /// state is equal to 0.
    fn len_gap(&self) -> PyResult<i32> {
        let length = self.coords.iter().filter(|x| **x < 0).collect::<Vec<&i32>>().len();
        Ok(length as i32)
    }
    
    // Format conversion

    #[staticmethod]
    /// Returns a linear space created using the given list of blocks.
    fn from_blocks(blocks: Vec<&Block>) -> PyResult<CoordSpace> {
        if blocks.len() == 0 {
            let coords: Vec<i32> = Vec::new();
            return Ok(CoordSpace{ coords })
        }
        match blocks_to_arrays(blocks) {
            Ok((data, ids)) => {
                let mut new_data: Vec<i32> = Vec::new();
                for i in 0..data.len() {
                    let x = data[i];
                    let id = &ids[i];
                    if id == "s" {
                        new_data.push(x);
                    } else {
                        new_data.push(-1);
                    }
                }
                Ok(CoordSpace { coords: new_data })
            },
            Err(x) => return Err(x)
        }
        
    }

    #[staticmethod]
    fn from_arrays(data: Vec<i32>, ids: Vec<String>) -> PyResult<CoordSpace> {
        if data.len() != ids.len() {
            return Err(exceptions::ValueError::py_err("lengths of data and ids do not match"))
        }
        if data.len() == 0 {
            let coords: Vec<i32> = Vec::new();
            return Ok(CoordSpace{ coords })
        }
        let mut coords: Vec<i32> = Vec::new();
        for i in 0..data.len() {
            let x = data[i];
            let id = &ids[i];
            if id == "s" {
                coords.push(x);
            } else if id == "g" {
                coords.push(-1);
            } else {
                return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", x)))
            }
        }
        Ok(CoordSpace{ coords })
    }

    /// Returns the linear space as a list of blocks.
    fn to_blocks(&self) -> PyResult<Vec<Block>> {
        if self.coords.len() == 0 {
            return Ok(Vec::new())
        }
        // Declare variables
        let mut blocks: Vec<Block> = Vec::new();
        let mut last_start: i32 = self.coords[0];
        let mut last_id: String = match self.coords[0] {
            x if x > 0 => "s".to_string(),
            x if x == -1 => "g".to_string(),
            x => return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", x))),
        };
        let mut negative_length: i32 = 0;

        for i in 1..self.coords.len() {
            let c_id: String = match self.coords[0] {
                x if x > 0 => "s".to_string(),
                x if x == -1 => "g".to_string(),
                x => return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", x))),
            };
            let c_pos = self.coords[i];
            let p_pos = self.coords[i-1];

            if c_pos == -1 && p_pos == -1 {
                negative_length += 1;
            } else if c_pos == -1 && p_pos > 0 {
                // Create new block and push
                blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
                // Assign current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
                negative_length = 0;
            } else if c_pos > 0 && p_pos == -1 {
                // Create new block and push
                blocks.push(Block{ id: last_id, start: 0, stop: negative_length});
                // Assign current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
                negative_length = 0;
            } else {
                if c_pos > p_pos + 1 {
                    // Create new block and push
                    blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
                    // Assgin current id as last_id and current pos as last_start
                    last_id = c_id;
                    last_start = c_pos;
                } else {
                    // Return an error
                    return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", c_pos)))
                }
            }
        }
        blocks.push(Block{ id: last_id, start: last_start, stop: self.coords.last().unwrap() + 1});
        Ok(blocks)
    }

    /// Returns the linear space as a list of integer coordinates.
    fn to_arrays(&self) -> PyResult<(Vec<i32>, Vec<String>)> {
        let coords = self.coords.clone();
        let mut ids: Vec<String> = Vec::new();
        for coord in self.coords.iter() {
            if *coord > 0 {
                ids.push("s".to_string());
            } else if *coord == -1 {
                ids.push("g".to_string())
            } else {
                return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", coord)))
            }
        }
        Ok((coords, ids))
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
    fn copy(&self) -> PyResult<CoordSpace> {
        let coords = self.coords.clone();
        Ok(CoordSpace{ coords })
    }

}

#[pyproto]
impl PyObjectProtocol for CoordSpace {
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
        Ok(format!("CoordSpace(start={}, stop={}, length={})", start, stop, length))
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

impl CoordSpace {

}

#[pyfunction]
/// blocks_to_arrays(block_list)
/// 
/// Converts a list of Blocks into an explicit listing of positions.
/// Returns a list of integers.
pub fn blocks_to_arrays(blocks: Vec<&Block>) -> PyResult<(Vec<i32>, Vec<String>)> {
    let mut data: Vec<i32> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    for Block{ id, start, stop} in blocks.iter() {
        let mut new_coords: Vec<i32> = (*start..*stop).collect();
        let mut new_ids: Vec<String> = vec![format!("{}", id); new_coords.len()];
        data.append(&mut new_coords);
        ids.append(&mut new_ids);
    }
    Ok((data, ids))
}

#[pyfunction]
/// arrays_to_blocks(data, ids)
/// 
/// Converts an explicit list of positions into a list of blocks.
/// Returns a list of Block objects.
pub fn arrays_to_blocks(data: Vec<i32>, ids: Vec<String>) -> PyResult<Vec<Block>> {
    if data.len() != ids.len() {
        return Err(exceptions::ValueError::py_err("lengths of data and ids do not match"))
    }
    if data.len() == 0 {
        return Ok(Vec::new())
    }
    // Declare variables
    let mut blocks: Vec<Block> = Vec::new();
    let mut last_id: &str = &ids[0];
    let mut last_start: i32 = data[0];

    for i in 1..data.len() {
        let c_id = &ids[i];
        let c_pos = data[i];
        let p_pos = data[i-1];
        // Scenarios:
        // 1a) current and last ids are the same, current +1 of previous
        // 1b) current and last ids are the same, current NOT +1 of previous
        // 2a) current and last ids NOT the same, current +1 of previous
        // 2b) current and last ids NOT the same, current NOT +1 of previous
        //
        // 2a and 2b are the same scenario, because change in ID should always
        // generate a new block
        if c_id == last_id {
            if c_pos != p_pos + 1 {
                // Create new block and push
                blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
                // Assgin current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
            }
        } else {
            // Create new block and push
            blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
            // Assign current id as last_id and current pos as last_start
            last_id = c_id;
            last_start = c_pos;
        }
    }
    blocks.push(Block{ id: last_id, start: last_start, stop: data.last().unwrap() + 1});
    Ok(blocks)
}

#[pymodinit]
fn position(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Block>()?;
    m.add_class::<CoordSpace>()?;

    Ok(())
}