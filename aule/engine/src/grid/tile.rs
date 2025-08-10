//! Tiling of the geodesic grid into compact tiles with halo rings.

use smallvec::SmallVec;

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
    /// Build a tiling for the provided grid using rectangular blocks per face and a halo of `halo_depth`.
    /// Construct a tiling. Current implementation is a trivial single-tile fallback.
    pub fn new(grid: &Grid, halo_depth: u8, target_tile_cells: usize) -> Self {
        // Fallback trivial tiling: one tile owns everything, halo empty. Real face-block tiler in follow-ups.
        let mut interior: Vec<u32> = (0..grid.cells as u32).collect();
        interior.sort_unstable();
        let tiles = vec![Tile {
            id: 0,
            interior: interior.clone(),
            halo: Vec::new(),
            neighbors: SmallVec::new(),
        }];
        let owner = vec![0u32; grid.cells];
        let local_index = (0..grid.cells as u32).collect();
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
