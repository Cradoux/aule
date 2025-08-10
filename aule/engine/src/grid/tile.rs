//! Tiling of the geodesic grid into compact tiles with halo rings.

use smallvec::SmallVec;
use std::collections::{BTreeSet, HashMap};

use super::Grid;

/// A single tile over the global grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Tile {
    /// Tile identifier
    pub id: u32,
    /// Interior cells owned by this tile
    pub interior: Vec<u32>,
    /// Halo cells (do not include interior)
    pub halo: Vec<u32>,
    /// Neighbor tile ids (interior adjacency), sorted unique
    pub neighbors: SmallVec<[u32; 6]>,
}

/// The full tiling description.
#[derive(Debug, Clone, PartialEq)]
pub struct Tiling {
    /// All tiles
    pub tiles: Vec<Tile>,
    /// Configured halo depth (k-rings captured)
    pub halo_depth: u8,
    /// Target interior size per tile
    pub target_tile_cells: usize,
    /// For each global cell id, the owning tile id (interior only)
    pub owner: Vec<u32>,
    /// For each global cell id, the local index within its owning tile interior
    pub local_index: Vec<u32>,
}

impl Tiling {
    /// Build a tiling for the provided grid using deterministic rectangular blocks per icosahedron face
    /// and a halo of `halo_depth` rings constructed via BFS on `n1` adjacency.
    pub fn new(grid: &Grid, halo_depth: u8, target_tile_cells: usize) -> Self {
        let f = grid.frequency.max(1);
        // Reconstruct per-face lattice -> global id mapping exactly as in Grid::new
        let per_face_grid_ids = reconstruct_per_face_grid_ids(f);

        // Determine block size in lattice coordinates aiming for ~target_tile_cells per tile.
        // Use square-ish blocks of side ~ sqrt(target), clamped to lattice extents.
        let side = ((target_tile_cells as f64).sqrt() as usize).max(1);

        let mut tiles: Vec<Tile> = Vec::new();
        let mut assigned: Vec<bool> = vec![false; grid.cells];

        for grid_ids in per_face_grid_ids.iter() {
            let f_usize = f as usize;
            let mut i0 = 0usize;
            while i0 <= f_usize {
                let mut j0 = 0usize;
                while j0 <= f_usize.saturating_sub(i0) {
                    // Build interior cells for this block by scanning the clipped rectangle
                    let mut interior: Vec<u32> = Vec::new();
                    let i_end = (i0 + side).min(f_usize);
                    for (i_idx, row) in grid_ids.iter().enumerate().skip(i0).take(i_end - i0 + 1) {
                        let max_j = f_usize - i_idx;
                        if j0 > max_j {
                            continue;
                        }
                        let j_end = (j0 + side).min(max_j);
                        for &gid in row.iter().skip(j0).take(j_end - j0 + 1) {
                            if !assigned[gid as usize] {
                                interior.push(gid);
                                assigned[gid as usize] = true;
                            }
                        }
                    }

                    if !interior.is_empty() {
                        interior.sort_unstable();
                        let tile_id = tiles.len() as u32;
                        // Halo via BFS up to halo_depth rings
                        let halo = build_halo(&interior, grid, halo_depth);
                        tiles.push(Tile {
                            id: tile_id,
                            interior,
                            halo,
                            neighbors: SmallVec::new(),
                        });
                    }

                    if j0 == f_usize.saturating_sub(i0) {
                        break;
                    }
                    j0 = j0.saturating_add(side + 1);
                }
                if i0 == f_usize {
                    break;
                }
                i0 = i0.saturating_add(side + 1);
            }
        }

        // Fallback in rare case some cells remain unassigned due to block stepping edges
        for gid in 0..grid.cells as u32 {
            if !assigned[gid as usize] {
                let tile_id = tiles.len() as u32;
                let interior = vec![gid];
                assigned[gid as usize] = true;
                let halo = build_halo(&interior, grid, halo_depth);
                tiles.push(Tile { id: tile_id, interior, halo, neighbors: SmallVec::new() });
            }
        }

        // Build owner/local_index from finalized interiors
        let mut owner: Vec<u32> = vec![u32::MAX; grid.cells];
        let mut local_index: Vec<u32> = vec![u32::MAX; grid.cells];
        for tile in &tiles {
            for (li, &gid) in tile.interior.iter().enumerate() {
                owner[gid as usize] = tile.id;
                local_index[gid as usize] = li as u32;
            }
        }

        // Build neighbor relations: if any interior cell of A has 1-ring neighbor owned by B ≠ A
        for tile in &mut tiles {
            let mut neigh: BTreeSet<u32> = BTreeSet::new();
            for &u in &tile.interior {
                for &v in &grid.n1[u as usize] {
                    let b = owner[v as usize];
                    if b != u32::MAX && b != tile.id {
                        neigh.insert(b);
                    }
                }
            }
            let mut sv = SmallVec::<[u32; 6]>::new();
            for t in neigh {
                sv.push(t);
            }
            tile.neighbors = sv;
        }

        Self { tiles, halo_depth, target_tile_cells, owner, local_index }
    }

    /// Return a view over a tile's interior and halo.
    /// Get a borrowed view for a specific tile id
    pub fn view(&self, tile_id: u32) -> TileView<'_> {
        let t = &self.tiles[tile_id as usize];
        TileView { interior: &t.interior, halo: &t.halo }
    }
}

/// Borrowed view of a tile's indices.
/// Borrowed view over a tile's interior and halo index slices.
pub struct TileView<'a> {
    /// Interior cell indices
    pub interior: &'a [u32],
    /// Halo cell indices
    pub halo: &'a [u32],
}

// ---- helpers ----

#[derive(Hash, Eq, PartialEq, Copy, Clone)]
struct EdgeKey(u32, u32, u32);

fn edge_key_from(f: u32, a_id: usize, b_id: usize, t_from_a: u32) -> (EdgeKey, u32) {
    let (lo_usize, hi_usize, t_from_lo) =
        if a_id < b_id { (a_id, b_id, t_from_a) } else { (b_id, a_id, f - t_from_a) };
    (EdgeKey(lo_usize as u32, hi_usize as u32, t_from_lo), t_from_lo)
}

fn reconstruct_per_face_grid_ids(f: u32) -> Vec<Vec<Vec<u32>>> {
    let (v0, faces) = icosahedron();
    // corner ids 0..11
    let mut corner_map: HashMap<usize, u32> = HashMap::new();
    for i in 0..v0.len() {
        corner_map.insert(i, i as u32);
    }
    let mut edge_map: HashMap<EdgeKey, u32> = HashMap::new();
    let mut next_id: u32 = v0.len() as u32;

    let mut per_face: Vec<Vec<Vec<u32>>> = Vec::with_capacity(20);
    for (fid, face) in faces.iter().enumerate() {
        let (a_id, b_id, c_id) = (face[0], face[1], face[2]);
        let mut grid_ids: Vec<Vec<u32>> = Vec::with_capacity(f as usize + 1);
        for i in 0..=f {
            let mut row: Vec<u32> = Vec::with_capacity((f - i) as usize + 1);
            for j in 0..=(f - i) {
                let k = f - i - j;
                let id = if j == 0 && k == 0 {
                    corner_map[&a_id]
                } else if i == 0 && k == 0 {
                    corner_map[&b_id]
                } else if i == 0 && j == 0 {
                    corner_map[&c_id]
                } else if k == 0 {
                    // AB edge, t=j
                    let t = j;
                    let (ekey, _) = edge_key_from(f, a_id, b_id, t);
                    if let Some(&eid) = edge_map.get(&ekey) {
                        eid
                    } else {
                        let eid = next_id;
                        edge_map.insert(ekey, eid);
                        next_id += 1;
                        eid
                    }
                } else if j == 0 {
                    // CA edge, t=i
                    let t = i;
                    let (ekey, _) = edge_key_from(f, c_id, a_id, t);
                    if let Some(&eid) = edge_map.get(&ekey) {
                        eid
                    } else {
                        let eid = next_id;
                        edge_map.insert(ekey, eid);
                        next_id += 1;
                        eid
                    }
                } else if i == 0 {
                    // BC edge, t=k
                    let t = k;
                    let (ekey, _) = edge_key_from(f, b_id, c_id, t);
                    if let Some(&eid) = edge_map.get(&ekey) {
                        eid
                    } else {
                        let eid = next_id;
                        edge_map.insert(ekey, eid);
                        next_id += 1;
                        eid
                    }
                } else {
                    // Interior unique to this face, assign fresh id in this pass order
                    let nid = next_id;
                    next_id += 1;
                    nid
                };
                row.push(id);
            }
            grid_ids.push(row);
        }
        let _ = fid; // avoid unused in some builds
        per_face.push(grid_ids);
    }
    per_face
}

fn build_halo(interior: &[u32], grid: &Grid, halo_depth: u8) -> Vec<u32> {
    if halo_depth == 0 {
        return Vec::new();
    }
    let mut visited: Vec<bool> = vec![false; grid.cells];
    for &u in interior {
        visited[u as usize] = true;
    }
    let mut frontier: Vec<u32> = interior.to_vec();
    let mut halo_set: BTreeSet<u32> = BTreeSet::new();
    for _ in 0..halo_depth {
        let mut next: Vec<u32> = Vec::new();
        for &u in &frontier {
            for &v in &grid.n1[u as usize] {
                let vi = v as usize;
                if !visited[vi] {
                    visited[vi] = true;
                    halo_set.insert(v);
                    next.push(v);
                }
            }
        }
        frontier = next;
    }
    // Ensure halo excludes interior (already ensured by visited flag) and sorted
    halo_set.into_iter().collect()
}

// Minimal icosahedron vertices/faces used to reconstruct lattice mapping
fn icosahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let phi = (1.0 + 5.0_f64.sqrt()) * 0.5;
    let a = 1.0;
    let b = 1.0 / phi;
    let normalize3 = |v: [f64; 3]| -> [f64; 3] {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / n, v[1] / n, v[2] / n]
    };
    let verts = vec![
        normalize3([-a, b, 0.0]),
        normalize3([a, b, 0.0]),
        normalize3([-a, -b, 0.0]),
        normalize3([a, -b, 0.0]),
        normalize3([0.0, -a, b]),
        normalize3([0.0, a, b]),
        normalize3([0.0, -a, -b]),
        normalize3([0.0, a, -b]),
        normalize3([b, 0.0, -a]),
        normalize3([b, 0.0, a]),
        normalize3([-b, 0.0, -a]),
        normalize3([-b, 0.0, a]),
    ];
    let faces = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    (verts, faces)
}
