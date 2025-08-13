//! T-455: High-F tiled mode plan and validators.
//!
//! Thin wrapper around `grid::tile::Tiling` to expose a stable API for per-tile
//! passes with halo rings and simple validation utilities.

use crate::grid::{self, Grid};

/// A single tile (interior + halo) for staging per-tile passes.
#[derive(Debug, Clone, PartialEq)]
pub struct Tile {
    /// Interior cell ids owned by this tile
    pub interior: Vec<u32>,
    /// Halo cell ids required by the stencil (read-only)
    pub halo: Vec<u32>,
}

/// Complete tiling plan metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct TilePlan {
    /// Tiles in deterministic order
    pub tiles: Vec<Tile>,
    /// Halo rings captured (1 for 1-ring stencils)
    pub halo_rings: u32,
    /// Safety cap for interior+halo cells per tile
    pub max_tile_cells: usize,
}

/// Builder parameters for tiling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TileBuildParams {
    /// Target interior size per tile (cells)
    pub target_cells_per_tile: usize,
    /// Halo rings depth (1 or 2)
    pub halo_rings: u32,
    /// Hard cap on total cells (interior+halo) per tile
    pub max_tile_cells: usize,
}

/// Errors from building/validating tile plans.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum TileError {
    /// The max_tile_cells cap is smaller than a tile's (interior+halo) requirement.
    #[error("max_tile_cells too small for smallest tile (need={0}, max={1})")]
    TooSmall(usize, usize),
    /// The halo did not include a required 1-ring neighbor for an interior cell.
    #[error("halo incomplete for cell {0} in tile {1}")]
    HaloIncomplete(u32, usize),
}

/// Build a deterministic `TilePlan` using face-lattice blocks and BFS halos.
pub fn build_tileplan(grid: &Grid, p: TileBuildParams) -> Result<TilePlan, TileError> {
    let halo_depth = (p.halo_rings as u8).max(0);
    let tiling = grid::tile::Tiling::new(grid, halo_depth, p.target_cells_per_tile);
    let mut tiles: Vec<Tile> = Vec::with_capacity(tiling.tiles.len());
    for t in &tiling.tiles {
        let total = t.interior.len() + t.halo.len();
        if total > p.max_tile_cells {
            return Err(TileError::TooSmall(total, p.max_tile_cells));
        }
        tiles.push(Tile { interior: t.interior.clone(), halo: t.halo.clone() });
    }
    Ok(TilePlan { tiles, halo_rings: p.halo_rings, max_tile_cells: p.max_tile_cells })
}

/// Validate that each interior cell's 1-ring neighborhood is covered by interior∪halo.
pub fn validate_tileplan(grid: &Grid, plan: &TilePlan, halo_rings: u32) -> Result<(), TileError> {
    let want = halo_rings.max(1);
    if plan.halo_rings < want {
        // Not an error type by itself; validation below will catch missing neighbors.
    }
    for (tid, t) in plan.tiles.iter().enumerate() {
        // Fast membership check
        let mut in_set = vec![false; grid.cells];
        for &u in &t.interior {
            in_set[u as usize] = true;
        }
        for &u in &t.halo {
            in_set[u as usize] = true;
        }
        for &u in &t.interior {
            for &v in &grid.n1[u as usize] {
                if !in_set[v as usize] {
                    return Err(TileError::HaloIncomplete(u, tid));
                }
            }
        }
    }
    Ok(())
}
