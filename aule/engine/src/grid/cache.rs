//! Binary cache read/write for `Grid`.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use super::Grid;

/// Errors related to grid cache I/O.
#[derive(thiserror::Error, Debug)]
pub enum GridError {
    /// Wrapper for standard I/O errors
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// Cache file has wrong magic or version
    #[error("bad magic or version")]
    BadHeader,
    /// Unexpected data length
    #[error("unexpected data length")]
    BadLength,
}
const MAGIC: &[u8; 9] = b"AULEGRID\0";
const VERSION: u32 = 1;

impl Grid {
    /// Save grid to a compact little-endian binary cache.
    pub fn save_cache<P: AsRef<Path>>(&self, path: P) -> Result<(), GridError> {
        let mut f = File::create(path)?;
        f.write_all(MAGIC)?;
        f.write_all(&VERSION.to_le_bytes())?;
        f.write_all(&(self.frequency).to_le_bytes())?;
        let cells_u32 = self.cells as u32;
        f.write_all(&cells_u32.to_le_bytes())?;

        // pos_xyz
        for p in &self.pos_xyz {
            for v in p {
                f.write_all(&v.to_le_bytes())?;
            }
        }
        // latlon
        for a in &self.latlon {
            for v in a {
                f.write_all(&v.to_le_bytes())?;
            }
        }
        // area
        for a in &self.area {
            f.write_all(&a.to_le_bytes())?;
        }

        // n1 with u8 length prefix per cell
        for neigh in &self.n1 {
            let len = u8::try_from(neigh.len()).unwrap_or(255);
            f.write_all(&[len])?;
            for &idx in neigh {
                f.write_all(&idx.to_le_bytes())?;
            }
        }
        // n2
        for neigh in &self.n2 {
            let len = u8::try_from(neigh.len()).unwrap_or(255);
            f.write_all(&[len])?;
            for &idx in neigh {
                f.write_all(&idx.to_le_bytes())?;
            }
        }

        Ok(())
    }

    /// Load grid from a binary cache file.
    pub fn load_cache<P: AsRef<Path>>(path: P) -> Result<Self, GridError> {
        let mut f = File::open(path)?;
        let mut magic = [0u8; 9];
        f.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(GridError::BadHeader);
        }
        let v = read_u32(&mut f)?;
        if v != VERSION {
            return Err(GridError::BadHeader);
        }
        let frequency = read_u32(&mut f)?;
        let cells = read_u32(&mut f)? as usize;

        let mut pos_xyz = vec![[0f32; 3]; cells];
        for p in &mut pos_xyz {
            for v in p.iter_mut() {
                *v = read_f32(&mut f)?;
            }
        }
        let mut latlon = vec![[0f32; 2]; cells];
        for a in &mut latlon {
            for v in a.iter_mut() {
                *v = read_f32(&mut f)?;
            }
        }
        let mut area = vec![0f32; cells];
        for a in &mut area {
            *a = read_f32(&mut f)?;
        }

        let mut n1 = Vec::with_capacity(cells);
        for _ in 0..cells {
            let mut len = [0u8; 1];
            f.read_exact(&mut len)?;
            let len = len[0] as usize;
            let mut v = smallvec::SmallVec::<[u32; 6]>::new();
            for _ in 0..len {
                v.push(read_u32(&mut f)?);
            }
            n1.push(v);
        }
        let mut n2 = Vec::with_capacity(cells);
        for _ in 0..cells {
            let mut len = [0u8; 1];
            f.read_exact(&mut len)?;
            let len = len[0] as usize;
            let mut v = smallvec::SmallVec::<[u32; 12]>::new();
            for _ in 0..len {
                v.push(read_u32(&mut f)?);
            }
            n2.push(v);
        }

        let mut grid = Self {
            cells,
            pos_xyz,
            latlon,
            area,
            n1,
            n2,
            frequency,
            east_hat: Vec::new(),
            north_hat: Vec::new(),
            lengths_n1_rad: Vec::new(),
            face_offsets: Vec::new(),
            face_vert_ids: Vec::new(),
            face_corners: Vec::new(),
        };
        // Populate precomputed bases/lengths after load
        grid.precompute_local_bases_and_lengths();
        Ok(grid)
    }
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32, GridError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn read_f32<R: Read>(r: &mut R) -> Result<f32, GridError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(f32::from_le_bytes(b))
}
