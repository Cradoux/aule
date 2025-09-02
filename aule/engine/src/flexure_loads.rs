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
    // Optionally subtract simple mean here is not possible without area; keep as-is.
    f
}

/// Simple effective sediment density as a function of sediment thickness (m).
/// Models porosity loss with burial: ρ_sed(thk) = ρ_min + (ρ_max − ρ_min) (1 − e^(−thk/H)).
/// Defaults: ρ_min=1600, ρ_max=2400 kg/m³, H=2000 m.
pub fn sediment_rho_eff(thk_m: f32) -> f32 {
    let rho_min = 1600.0f32;
    let rho_max = 2400.0f32;
    let h = 2000.0f32;
    let t = thk_m.max(0.0);
    rho_min + (rho_max - rho_min) * (1.0 - (-t / h).exp())
}

/// Assemble load including explicit sediment contribution with compaction.
///
/// Components per cell i (with z = η − depth):
/// - Ocean water: ρ_w g (h_ocean − sed_submerged)
/// - Emergent crust (excl. subaerial sediments): ρ_c g max(z − sed_subaerial, 0)
/// - Subaerial sediments: ρ_sed_eff g sed_subaerial
/// - Submerged sediments: (ρ_sed_eff − ρ_w) g sed_submerged
pub fn assemble_load_with_sediments(
    grid: &Grid,
    depth_m: &[f32],
    sediment_m: &[f32],
    p: &LoadParams,
) -> Vec<f32> {
    assert_eq!(depth_m.len(), grid.cells);
    let mut f = vec![0.0f32; grid.cells];
    let rho_w_g = p.rho_w * p.g;
    let rho_c_g = p.rho_c * p.g;
    for i in 0..grid.cells {
        let d = depth_m[i];
        let z = p.sea_level_m - d; // elevation
        let h_ocean = (d - p.sea_level_m).max(0.0);
        let sed_thk = sediment_m.get(i).copied().unwrap_or(0.0).max(0.0);
        // Partition sediment into submerged vs subaerial by available water/land thickness
        let sed_submerged = sed_thk.min(h_ocean);
        let sed_subaerial = (sed_thk - sed_submerged).max(0.0);
        // Emergent crustal elevation excluding subaerial sediment top layer
        let elev_no_sed = (z - sed_subaerial).max(0.0);
        let rho_sed_g = sediment_rho_eff(sed_thk) * p.g;
        // Water over remaining ocean (not occupied by submerged sediment)
        let ocean_water = rho_w_g * (h_ocean - sed_submerged).max(0.0);
        // Emergent crust
        let crust_emergent = rho_c_g * elev_no_sed;
        // Sediments: subaerial fully add; submerged add relative to water they displace
        let sed_air = rho_sed_g * sed_subaerial;
        let sed_below = (rho_sed_g - rho_w_g) * sed_submerged;
        f[i] = ocean_water + crust_emergent + sed_air + sed_below;
    }
    f
}
