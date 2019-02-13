use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};


#[pyclass(subclass)]
#[derive(Clone)]
/// Block(id, start, stop, /)
/// --
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

    /// in_block(i)
    /// 
    /// Returns True if a given position is inside the block.
    /// Otherwise, returns False.
    fn in_block(&self, i: i32) -> PyResult<bool> {
        if i >= self.start && i < self.stop {
            return Ok(true)
        }
        Ok(false)
    }

    /// to_array()
    ///
    /// Converts the block into a list of positions.
    fn to_array(&self) -> PyResult<Vec<i32>> {
        Ok(self._to_array())
    }

    // Formatting methods

    /// to_compressed_str()
    ///
    /// Converts block into a compressed string representation.
    fn to_compressed_str(&self) -> PyResult<String> {
        self.__str__()
    }

    /// to_extended_str()
    ///
    /// Converts block into an extended (human-readable) string
    /// representation.
    fn to_extended_str(&self) -> PyResult<String> {
        Ok(format!("{}={}:{}", self.id, self.start, self.stop))
    }
    
    // TODO: Add a method to convert to CIGAR string
    // fn to_cigar_str(&self) -> PyResult<String> {
    // }

    /// to_array_str()
    /// 
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
/// LinearSpace(start, stop, init_state, /)
/// --
/// 
/// LinearSpace represents a discrete linear space stored
/// in intervals called Blocks.
pub struct LinearSpace {

    coords: Vec<[i32; 3]>,

}

#[pyproto]
impl PyObjectProtocol for LinearSpace {
    fn __repr__(&self) -> PyResult<String> {
        let lb = match self.lb() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let ub = match self.ub() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        let length = match self.len() {
            Ok(x) => x,
            Err(x) => return Err(x)
        };
        Ok(format!("LinearSpace(lb={}, ub={}, length={})", lb, ub, length))
    }
    
    fn __str__(&self) -> PyResult<String> {
        self.to_extended_str()
    }
}

#[pymethods]
impl LinearSpace {
    #[new]
    /// Creates a new LinearSpace object from an init_state, and start and stop
    /// coordinates.
    fn __new__(obj: &PyRawObject, start: i32, stop: i32, init_state: i32) -> PyResult<()> {
        if start > stop {
            return Err(exceptions::ValueError::py_err(
                format!("start must be less than stop: {} !< {}",
                        start, stop)))
        }
        obj.init(|_| {
            LinearSpace {
                coords: vec![[start, stop, init_state]],
            }
        })
    }

    /// extract(positions, /)
    /// --
    /// 
    /// Extracts coordinates according to its current positions
    /// as a new LinearSpace.
    fn extract(&self, positions: Vec<i32>) -> PyResult<LinearSpace> {
        if let Some(max) = positions.iter().max() {
            if *max >= self.coords.len() as i32 {
                return Err(exceptions::IndexError::py_err(
                    format!("index out of range: {}", max)))
            }
            // Unroll blocks into a vector of i32
            let (coord_list, id_list) = self.to_arrays()?;
            // Extract
            let mut ext_coord_list: Vec<i32> = Vec::new();
            let mut ext_id_list: Vec<i32> = Vec::new();
            for i in positions.iter() {
                ext_coord_list.push(coord_list[*i as usize]);
                ext_id_list.push(id_list[*i as usize]);
            }
            // Reassemble to blocks
            arrays_to_linspace(ext_coord_list, ext_id_list)
        } else {
            Ok(LinearSpace{ coords: self.coords.clone() })
        }
    }

    // Realtive position methods

    /// remove(positions, /)
    /// --
    /// 
    /// Removes points in linear space given based on a list of
    /// current positions.
    fn remove(&mut self, positions: Vec<i32>) -> PyResult<()> {
        // Check if positions list is empty or not using max()
        if let Some(max) = positions.iter().max() {
            // If largest relative position is larger than the total length
            // of the linear space, then return an error
            let length = self.len().unwrap();
            if *max >= length {
                return Err(exceptions::ValueError::py_err(
                    format!("index out of range: {}", max)))
            }

            // Create a vector of arrays representing relative start and stops
            let mut start_offset = 0;
            let rel_blocks: Vec<[i32;2]> = self.coords.iter()
                .map(|[a, z, _]| {
                    let length = z - a;
                    start_offset += length;
                    [start_offset-length, start_offset]
            }).collect();

            // Sort positions in reverse order
            let mut positions = positions;
            positions.sort_unstable();
            positions.reverse();
            
            // Iterate over positions in reverse
            let mut offset = 0;
            for &rel_pos in positions.iter() {
                let length = rel_blocks.len() - offset;
                for i in (0..length).rev() {
                    let [rel_start, rel_stop]: [i32; 2] = rel_blocks[i];
                    let [abs_start, abs_stop, id] = self.coords[i];
                    if rel_pos > rel_start && rel_pos < rel_stop - 1 {
                        // Remove block currently at j and get values
                        // Split this block at pos
                        // 10,11,12,13,14,15 : remove at index 2 (3rd pos)
                        // [_,10,16]
                        // 10,11,   13,14,15
                        // [_,10,12], [13:16]
                        // Insert two new blocks at j and j+1
                        let abs_pos = abs_start + (rel_pos - rel_start);
                        self.coords[i] = [abs_start, abs_pos, id];
                        self.coords.insert(i+1, [abs_pos+1, abs_stop, id]);
                        offset += 1;
                        break;
                    } else if rel_pos == rel_start {
                        // target is the start of the block
                        if abs_start+1 == abs_stop {
                            let _ = self.coords.remove(i);
                            // no offset increment
                        } else {
                            self.coords[i] = [abs_start+1, abs_stop, id];
                            offset += 1;
                        }
                        break;
                    } else if rel_pos == rel_stop - 1 {
                        // target is the end of the block
                        self.coords[i] = [abs_start, abs_stop-1, id];
                        // no offset increment
                        break;
                    }
                }
            }
        } 
        Ok(())
    }

    /// retain(positions, /)
    /// 
    /// Retains points in linear space specified by a
    /// list of positions to keep.
    fn retain(&mut self, positions: Vec<i32>) -> PyResult<()> {
        // Check if positions list is empty or not using max()
        if let Some(max) = positions.iter().max() {
            // If largest relative position is larger than the total length
            // of the linear space, then return an error
            let length = self.len().unwrap();
            if *max >= length {
                return Err(exceptions::ValueError::py_err(
                    format!("index out of range: {}", max)))
            }
            let inverse_rel_positions: Vec<i32> = (0..length)
                .filter(|x| !positions.contains(x))
                .collect();
            return self.remove(inverse_rel_positions)
        } 
        Ok(())
    }

    // Absolute position methods

    // /// remove_abs(positions, /)
    // /// --
    // /// 
    // /// Removes points in linear space given based on a list of coordinates.
    // fn remove_abs(&mut self, positions: Vec<i32>) -> PyResult<()> {
    //     let mut positions = positions;
    //     positions.sort_unstable();
    //     positions.reverse();
    //     let mut offset = 0;
    //     for &pos in positions.iter() {
    //         let length = self.coords.len() - offset;
    //         for i in (0..length).rev() {
    //             let [start, stop, id]: [i32; 3] = self.coords[i];
    //             if pos > start && pos < stop - 1 {
    //                 // Remove block currently at j and get values
    //                 // Split this block at pos
    //                 // 10,11,12,13,14,15 : remove at index 2 (3rd pos)
    //                 // [_,10,16]
    //                 // 10,11,   13,14,15
    //                 // [_,10,12], [13:16]
    //                 // Insert two new blocks at j and j+1
    //                 self.coords[i] = [start, pos, id];
    //                 self.coords.insert(i+1, [pos+1, stop, id]);
    //                 offset += 1;
    //             } else if pos == start {
    //                 if start == stop - 1 {
    //                     let _ = self.coords.remove(i);
    //                     // no offset increment
    //                 } else {
    //                     self.coords[i] = [pos+1, stop, id];
    //                     // no offset increment
    //                 }
    //             } else if pos == stop - 1 {
    //                 self.coords[i] = [start, pos, id];
    //                 // no offset increment
    //             }
    //         }
    //     }
    //     Ok(())
    // }

    // /// retain_abs(coords, /)
    // /// --
    // /// 
    // /// Retains points in linear space specified by a
    // /// list of coordinates to keep.
    // fn retain_abs(&mut self, coords: Vec<i32>) -> PyResult<()> {
    //     if let Some([_, _, stop]) = self.coords.last() {
    //         let inverse_ilist: Vec<i32> = (0..*stop)
    //                                         .filter(|x| !coords.contains(x))
    //                                         .collect();
    //         return self.remove_abs(inverse_ilist)
    //     };
    //     Err(exceptions::ValueError::py_err("cannot perform retain on \
    //                                         dimension: block list is empty"))
    // }

    // TODO: Add insert and append methods

    // start, stop, full_len

    /// lb()
    /// --
    /// 
    /// Returns the lower bound of the linear space.
    fn lb(&self) -> PyResult<i32> {
        match self.coords.first() {
            Some([x, _, _]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err(
                "linear space is empty"))
        }
    }

    /// ub()
    /// --
    /// 
    /// Returns the upper bound of the linear space.
    /// This value is not part of the space.
    fn ub(&self) -> PyResult<i32> {
        match self.coords.last() {
            Some([_, x, _]) => Ok(*x),
            None => return Err(exceptions::ValueError::py_err(
                "linear space is empty"))
        }
    }

    /// Returns the total length of the linear space.
    fn len(&self) -> PyResult<i32> {
        if self.coords.len() == 0 {
            return Ok(0)
        }
        let mut length = 0;
        for [start, stop, _] in self.coords.iter() {
            length += stop - start;
        }
        Ok(length)
    }
    
    // Format conversion

    /// to_blocks()
    /// --
    /// 
    /// Returns the linear space as a list of blocks.
    fn to_blocks(&self) -> PyResult<Vec<Block>> {
        if self.coords.len() == 0 {
            return Ok(Vec::new())
        }
        let mut blocks: Vec<Block> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            blocks.push(Block{ id: format!("{}", id), start: *start, stop: *stop });
        }
        Ok(blocks)
    }

    // /// Returns the linear space as a list of point coordinates.
    // fn to_points(&self) -> PyResult<Vec<Point>> {
    // }

    /// to_list()
    /// --
    /// 
    /// Returns the linear space as a list of start, stop, and id tuples.
    fn to_list(&self) -> PyResult<Vec<(i32, i32, i32)>> {
        if self.coords.len() == 0 {
            return Ok(Vec::new())
        }
        let mut list: Vec<(i32, i32, i32)> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            list.push((*start, *stop, *id));
        }
        Ok(list)
    }

    /// to_arrays()
    /// --
    /// 
    /// Returns the linear space as corresponding coordinates and id lists.
    fn to_arrays(&self) -> PyResult<(Vec<i32>, Vec<i32>)> {
        if self.coords.len() == 0 {
            return Ok((Vec::new(), Vec::new()))
        }
        let mut coords: Vec<i32> = Vec::new();
        let mut ids: Vec<i32> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            for i in *start..*stop {
                coords.push(i);
                ids.push(*id);
            }
        }
        Ok((coords, ids))
    }

    // Formatting methods

    /// to_compressed_str()
    /// --
    /// 
    /// Converts block into a compressed string representation.
    fn to_compressed_str(&self) -> PyResult<String> {
        if self.coords.len() == 0 {
            return Ok(String::new())
        }
        let mut strings: Vec<String> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            strings.push(format!("{}={}", id, stop-start));
        }
        Ok(strings.join(";"))
    }

    /// to_extended_str()
    /// --
    /// 
    /// Converts block into an extended (human-readable) string
    /// representation.
    fn to_extended_str(&self) -> PyResult<String> {
        if self.coords.len() == 0 {
            return Ok(String::new())
        }
        let mut strings: Vec<String> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            strings.push(format!("{}={}:{}", id, start, stop));
        }
        Ok(strings.join(";"))
    }

    /// Converts block into comma-separated list of positions.
    fn to_array_str(&self) -> PyResult<String> {
        if self.coords.len() == 0 {
            return Ok(String::new())
        }
        let mut strings: Vec<String> = Vec::new();
        for [start, stop, id] in self.coords.iter() {
            let mut b_strings: Vec<String> = Vec::new();
            b_strings.push(format!("{}=", id));
            for i in *start..*stop {
                b_strings.push(format!("{}", i));
            }
            strings.push(b_strings.join(","));
        }
        Ok(strings.join(";"))  
    }

    /// copy()
    /// --
    /// 
    /// Returns a deep copy of the current linear space.
    fn copy(&self) -> PyResult<LinearSpace> {
        let coords = self.coords.clone();
        Ok(LinearSpace{ coords })
    }
}


#[pyfunction]
/// blocks_to_linspace(blocks, /)
/// --
/// 
/// Returns a linear space created using the given list of blocks.
pub fn blocks_to_linspace(blocks: Vec<&Block>) -> PyResult<LinearSpace> {
    let mut coords: Vec<[i32; 3]> = Vec::new();
    for Block{ id, start, stop } in blocks.iter() {
        if let Ok(v) = id.parse::<i32>() {
            coords.push([*start, *stop, v]);
        } else {
            return Err(exceptions::ValueError::py_err("error converting Block id to int"))
        }
    }
    Ok(LinearSpace{ coords })
}

#[pyfunction]
/// list_to_linspace(coords, /)
/// --
/// 
/// Returns a linear space created using the given coordinate list.
pub fn list_to_linspace(coords: Vec<(i32, i32, i32)>) -> PyResult<LinearSpace> {
    let coords: Vec<[i32; 3]> = coords.iter().map(|(a, b, c)| [*a, *b, *c]).collect();
    Ok(LinearSpace{ coords })
}

#[pyfunction]
/// arrays_to_linspace(coords, ids, /)
/// --
/// 
/// Returns a linear space based on the corresponding lists of coordinates
/// and ids.
pub fn arrays_to_linspace(coords: Vec<i32>, ids: Vec<i32>) -> PyResult<LinearSpace> {
    if coords.len() != ids.len() {
        return Err(exceptions::ValueError::py_err("lists do not have the same length"))
    }
    let mut new_coords: Vec<[i32; 3]> = Vec::new();
    let mut last_id = ids[0];
    let mut last_start = coords[0];
    for i in 1..ids.len() {
        let c_id = ids[i];
        let c_pos = coords[i];
        let p_pos = coords[i-1];
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
                new_coords.push([last_start, p_pos + 1, last_id]);
                // Assgin current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
            }
        } else {
            // Create new block and push
            new_coords.push([last_start, p_pos + 1, last_id]);
            // Assign current id as last_id and current pos as last_start
            last_id = c_id;
            last_start = c_pos;
        }
    }
    new_coords.push(
        [last_start, coords.last().unwrap() + 1, last_id]);
    Ok(LinearSpace{ coords: new_coords })
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

    /// extract(coordinates)
    /// 
    /// Extracts coordinates by relative positions as a new CoordSpace.
    fn extract(&self, coords: Vec<i32>) -> PyResult<CoordSpace> {
        if let Some(max) = coords.iter().max() {
            if *max >= self.coords.len() as i32 {
                return Err(exceptions::IndexError::py_err(format!("index out of range: {}", max)))
            }
            let mut new_coords: Vec<i32> = Vec::new();
            for i in coords.iter() {
                new_coords.push(self.coords[*i as usize]);
            }
            Ok(CoordSpace{ coords: new_coords })
        } else {
            Ok(CoordSpace { coords: self.coords.clone()})
        }
    }

    /// remove(coordinates)
    /// 
    /// Removes points in linear space given based on a list of relative
    /// coordinates.
    fn remove(&mut self, coords: Vec<i32>) -> PyResult<()> {
        if let Some(max) = coords.iter().max() {
            if *max >= self.coords.len() as i32 {
                return Err(exceptions::IndexError::py_err(format!("index out of range: {}", max)))
            }
            self.coords = self.coords.iter().enumerate().filter(|(i, _)| !coords.contains(&(*i as i32))).map(|(_, x)| *x ).collect();
            Ok(())
        } else {
            Ok(())
        }
        
    }

    /// retain(coordinates)
    /// 
    /// Retains points in linear space specified by a
    /// list of coordinates to keep.
    fn retain(&mut self, coords: Vec<i32>) -> PyResult<()> {
        if let Some(max) = coords.iter().max() {
            if *max >= self.coords.len() as i32 {
                return Err(exceptions::IndexError::py_err(format!("index out of range: {}", max)))
            }
            self.coords = self.coords.iter().enumerate().filter(|(i, _)| coords.contains(&(*i as i32))).map(|(_, x)| *x ).collect();
            Ok(())
        } else {
            self.coords = Vec::new();
            Ok(())
        }        
    }

    // /// Inserts into the linear space at the given position.
    // fn insert(&mut self, pos: i32, start: i32, length: i32) -> PyResult<()> {
    //     // Insert to start of list if pos is 0,
    //     // Append to end of list if pos is the length,
    //     // Otherwise split at the position and concat
    //     if pos == 0 {
    //         let mut new_coords: Vec<i32> = (start..start+length).collect();
    //         new_coords.append(&mut self.coords);
    //         self.coords = new_coords;
    //         return Ok(())
    //     } else if pos == self.coords.len() as i32 {
    //         return self.append(start, length)
    //     } // TODO: else if pos > self.coords.len()
    //     // Split list into two and combine
    //     let mut new_coords: Vec<i32> = (start..start+length).collect();
    //     let (left, right) = self.coords.split_at(pos as usize);
    //     let mut left: Vec<i32> = left.iter().map(|x| *x ).collect();
    //     let mut right: Vec<i32> = right.iter().map(|x| *x ).collect();
    //     left.append(&mut new_coords);
    //     left.append(&mut right);
    //     self.coords = left;
    //     Ok(())
    // }

    // /// Appends to the end of the linear space.
    // fn append(&mut self, start: i32, length: i32) -> PyResult<()> {
    //     let mut new_coords: Vec<i32> = (start..start+length).collect();
    //     self.coords.append(&mut new_coords);
    //     Ok(())
    // }

    // start, stop, full_len

    /// start()
    /// 
    /// Returns the value at the start of the linear space.
    fn start(&self) -> PyResult<i32> {
        Ok(self.coords[0])
    }

    /// stop()
    /// 
    /// Returns the end value of the linear space.
    /// This value is not part of the space.
    fn stop(&self) -> PyResult<i32> {
        if self.coords.len() > 0 {
            Ok(self.coords[self.coords.len() - 1])
        } else {
            Ok(0)
        }
    }

    /// len_all()
    /// 
    /// Returns the total length of the linear space.
    fn len_all(&self) -> PyResult<i32> {
        Ok(self.coords.len() as i32)
    }

    /// len_seq()
    /// 
    /// Returns the total length of the linear space where the
    /// state is equal to 1.
    fn len_seq(&self) -> PyResult<i32> {
        let length = self.coords.iter().filter(|x| **x > 0).collect::<Vec<&i32>>().len();
        Ok(length as i32)
    }

    /// len_gap()
    /// 
    /// Returns the total length of the linear space where the
    /// state is equal to 0.
    fn len_gap(&self) -> PyResult<i32> {
        let length = self.coords.iter().filter(|x| **x < 0).collect::<Vec<&i32>>().len();
        Ok(length as i32)
    }
    
    // Format conversion

    #[staticmethod]
    /// from_blocks(blocks)
    /// 
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
                    } else if id == "g" {
                        new_data.push(-1);
                    } else {
                        return Err(exceptions::ValueError::py_err(format!("unsupported ID: {}. Use \"s\" for sequence or \"g\" for gap.", id)))
                    }
                }
                Ok(CoordSpace { coords: new_data })
            },
            Err(x) => return Err(x)
        }
        
    }

    #[staticmethod]
    /// from_arrays(coordinates, ids)
    /// 
    /// Returns a linear space created using the corresponding lists of
    /// coordinates and ids.
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
                return Err(exceptions::ValueError::py_err(format!("unsupported ID: {}. Use \"s\" for sequence or \"g\" for gap.", id)))
            }
        }
        Ok(CoordSpace{ coords })
    }

    /// to_blocks()
    /// 
    /// Returns the linear space as a list of blocks.
    fn to_blocks(&self) -> PyResult<Vec<Block>> {
        if self.coords.len() == 0 {
            return Ok(Vec::new())
        }
        // Declare variables
        let mut blocks: Vec<Block> = Vec::new();
        let mut last_start: i32 = self.coords[0];
        let mut last_id: String = match self.coords[0] {
            x if x >= 0 => "s".to_string(),
            x if x == -1 => "g".to_string(),
            x => return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", x))),
        };
        let mut negative_length: i32 = 0;

        for i in 1..self.coords.len() {
            let c_id: String = match self.coords[0] {
                x if x >= 0 => "s".to_string(),
                x if x == -1 => "g".to_string(),
                x => return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", x))),
            };
            let c_pos = self.coords[i];
            let p_pos = self.coords[i-1];

            if c_pos == -1 && p_pos == -1 {
                negative_length += 1;
            } else if c_pos < -1 || p_pos < -1 {
                // Return an error
                return Err(exceptions::ValueError::py_err(format!("unexpected coordinate value: {}", c_pos)))
            } else if c_pos == -1 && p_pos >= 0 {
                // Create new block and push
                blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
                // Assign current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
                negative_length = 0;
            } else if c_pos >= 0 && p_pos == -1 {
                // Create new block and push
                blocks.push(Block{ id: last_id, start: 0, stop: negative_length});
                // Assign current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
                negative_length = 0;
            } else if c_pos >= 0 && p_pos >= 0 {
                if c_pos != p_pos + 1 {
                    // Create new block and push
                    blocks.push(Block{ id: last_id, start: last_start, stop: p_pos + 1});
                    // Assgin current id as last_id and current pos as last_start
                    last_id = c_id;
                    last_start = c_pos;
                }
            }
        }
        blocks.push(Block{ id: last_id, start: last_start, stop: self.coords.last().unwrap() + 1});
        Ok(blocks)
    }

    /// to_arrays()
    /// 
    /// Returns the linear space as a list of integer coordinates.
    fn to_arrays(&self) -> PyResult<(Vec<i32>, Vec<String>)> {
        let coords = self.coords.clone();
        let mut ids: Vec<String> = Vec::new();
        for coord in self.coords.iter() {
            if *coord >= 0 {
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

    /// to_compressed_str()
    /// 
    /// Converts block into a compressed string representation.
    fn to_compressed_str(&self) -> PyResult<String> {
        self.__str__()
    }

    /// to_extended_str()
    /// 
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

    /// to_array_str()
    /// 
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

    /// copy()
    /// 
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
        match self.to_blocks() {
            Ok(blocks) => {
                for block in blocks {
                    if let Ok(s) = block.__str__() {
                        strings.push(s);
                    } else {
                        return Err(exceptions::ValueError::py_err("cannot get string representation of block"))
                    }
                }
            },
            Err(_) => {
                return Err(exceptions::ValueError::py_err(format!("cannot generate blocks")))
            }
            
        }
        Ok(strings.join(","))
    }
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
                blocks.push(Block{ id: last_id.to_string(), start: last_start, stop: p_pos + 1});
                // Assgin current id as last_id and current pos as last_start
                last_id = c_id;
                last_start = c_pos;
            }
        } else {
            // Create new block and push
            blocks.push(Block{ id: last_id.to_string(), start: last_start, stop: p_pos + 1});
            // Assign current id as last_id and current pos as last_start
            last_id = c_id;
            last_start = c_pos;
        }
    }
    blocks.push(Block{ id: last_id.to_string(), start: last_start, stop: data.last().unwrap() + 1});
    Ok(blocks)
}

#[pymodinit]
fn position(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Block>()?;
    m.add_class::<LinearSpace>()?;
    m.add_class::<CoordSpace>()?;

    m.add_function(wrap_function!(blocks_to_linspace))?;
    m.add_function(wrap_function!(list_to_linspace))?;
    m.add_function(wrap_function!(arrays_to_linspace))?;

    Ok(())
}
