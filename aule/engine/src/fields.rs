//! Field views and tile-aware accessors (index indirection only).

use crate::grid::tile::TileView;

/// Example of a structure-of-arrays storage for a scalar field.
#[derive(Debug)]
pub struct ScalarField {
    /// Global SoA storage (length equals number of grid cells)
    pub data: Vec<f32>,
}

impl ScalarField {
    /// Create a new scalar field with `len` cells initialized to 0.0
    pub fn new(len: usize) -> Self {
        Self { data: vec![0.0; len] }
    }

    /// Get a tile view as slices into the global storage by indices.
    /// Build a tile-aware view that carries index slices into the global storage.
    pub fn view<'a>(&'a self, tv: TileView<'a>) -> ScalarTileView<'a> {
        ScalarTileView { interior: tv.interior, halo: tv.halo, data: &self.data }
    }
}

/// A tile view for a scalar field that holds index slices and a reference to the data.
/// A borrowed tile view for a scalar field holding index slices and data reference.
pub struct ScalarTileView<'a> {
    /// Interior cell indices for this tile
    pub interior: &'a [u32],
    /// Halo cell indices (excludes interior)
    pub halo: &'a [u32],
    /// Reference to the global scalar data buffer
    pub data: &'a [f32],
}
