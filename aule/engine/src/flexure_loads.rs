//! Flexure load assembly from depth/bathymetry.
//!
//! Given per-cell depths in meters (positive downward), we form a vertical surface load `f` (N/m^2):
//!
//! [ f = \rho_w g \max(\text{depth}, 0) + \rho_c g \max(-\text{depth}, 0) ]
//!
//! where water exerts pressure over oceans (depth > 0) and continental crust weight applies over land
//! (depth ≤ 0, with elevation = -depth). This keeps continuity across sea level.

use crate::grid::Grid;

/// Inputs in SI. `depth_m > 0`: ocean depth; `depth_m ≤ 0`: land elevation = `-depth_m`.
#[derive(Clone, Copy)]
pub struct LoadParams {
    /// Water density (kg/m^3), e.g. 1030
    pub rho_w: f32,
    /// Continental crust density (kg/m^3), e.g. 2900
    pub rho_c: f32,
    /// Gravitational acceleration (m/s^2), e.g. 9.81
    pub g: f32,
    /// Sea level reference (meters). Typically 0; if a solved offset is available, pass it here.
    pub sea_level_m: f32,
}

/// Assemble vertical surface load (N/m^2) per cell from depth/elevation.
///
/// Formula: `f = rho_w * g * max(depth, 0) + rho_c * g * max(-depth, 0)` where `depth` is measured
/// positive downward relative to `sea_level_m`.
pub fn assemble_load_from_depth(grid: &Grid, depth_m: &[f32], p: &LoadParams) -> Vec<f32> {
    assert_eq!(depth_m.len(), grid.cells);
    let mut f = vec![0.0f32; grid.cells];
    let rho_w_g = p.rho_w * p.g;
    let rho_c_g = p.rho_c * p.g;
    for i in 0..grid.cells {
        let d_rel = depth_m[i] - p.sea_level_m;
        // ocean pressure for depth>0, crustal load for elevation>0 (depth≤0)
        let ocean = if d_rel > 0.0 { d_rel } else { 0.0 };
        let land = if d_rel < 0.0 { -d_rel } else { 0.0 };
        f[i] = rho_w_g * ocean + rho_c_g * land;
    }
    f
}
