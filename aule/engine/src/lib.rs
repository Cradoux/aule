//! Aulë engine crate stub.
//! Minimal; includes GPU helpers for tests.
#![deny(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

/// Steady-state oceanic age and bathymetry mapping.
pub mod age;
/// Plate boundary classification.
pub mod boundaries;
pub use boundaries::{Boundaries, BoundaryStats, EdgeClass, EdgeKin};
/// CFL (Courant-Friedrichs-Lewy) limiter for numerical stability.
pub mod cfl;
/// Configuration types shared between pipeline and world modules.
pub mod config;
/// Centralized elevation buffer management with caching.
pub mod elevation_manager;
/// Flexure load assembly helpers.
pub mod flexure_loads;

/// Arc/terrane accretion at O–C margins.
pub mod accretion;
/// Synthetic continental mask generator and applier.
pub mod continent;
/// Erosion helpers (stream-power law, slopes, drainage area).
pub mod erosion;
/// CPU face/tri picker using shared `aule-geo`.
pub mod face_pick;
/// Field and tiling views.
pub mod fields;
/// 1D elastic plate flexure (CPU reference; analytic + CG solver).
pub mod flexure;
/// GPU flexure multigrid scaffold (WGSL).
pub mod flexure_gpu;
/// Flexible flexure backend management with CPU/GPU selection.
pub mod flexure_manager;
/// Force-balance update for Euler poles.
pub mod force_balance;
/// Geometric utilities used by solvers.
pub mod geo;
/// Minimal GPU helper for tests (T-020).
pub mod gpu;
/// Geodesic grid module.
pub mod grid;
/// Hypsometry helpers.
pub mod hypsometry;
/// Global sea-level (isostasy MVP) utilities.
pub mod isostasy;
/// Orogeny (C–C collision) bilateral thickening/uplift.
pub mod orogeny;
/// Plate adjacency and triple-junction diagnostics.
pub mod plate_network;
/// Plates seeding and velocities.
pub mod plates;
/// Ridge births and fringe assignment (CPU pass).
pub mod ridge;
/// Continental rifting and passive margins.
pub mod rifting;
/// Long-term eustatic sea-level controller.
pub mod sea_level;
/// Sediment transport/deposition helpers.
pub mod sediment;
/// Semi-Lagrangian advection helpers (CPU, deterministic).
pub mod sl_advect;
/// Snapshot IO utilities (CSV/binary)
pub mod snapshots;
/// World stepper.
pub mod stepper;
/// Subduction bands and bathymetry adjustments (CPU pass).
pub mod subduction;
/// Fluvial erosion, hillslope diffusion, and sediment transport/deposition orchestrator.
pub mod surface;
/// High-F tiling plan (T-455)
pub mod tileplan;
pub use world::{SimplePreset, SimpleReport};
/// Centralized physics pipeline for Simple/Advanced.
pub mod pipeline;
/// Transform pull-apart/restraining bands (CPU pass).
pub mod transforms;
/// Units-of-measure lightweight wrappers.
pub mod units;
/// Unified pipeline that replaces both world::step_once and pipeline::step_full.
pub mod unified_pipeline;
/// Small utilities.
pub mod util;
/// World state.
pub mod world;

/// Centralized physical constants (SI units) for densities, gravity, and process limits.
#[derive(Clone, Copy, Debug)]
pub struct PhysConsts {
    /// Water density (kg/m^3)
    pub rho_w_kg_per_m3: f32,
    /// Continental crust density (kg/m^3)
    pub rho_c_kg_per_m3: f32,
    /// Mantle density (kg/m^3)
    pub rho_m_kg_per_m3: f32,
    /// Air density (kg/m^3)
    pub rho_air_kg_per_m3: f32,
    /// Gravity (m/s^2)
    pub g_m_per_s2: f32,
    /// Reference continental crust thickness (m)
    pub th_ref_continental_m: f32,
    /// Earth radius (m)
    pub r_earth_m: f64,
    /// Maximum erosion rate per Myr (m/Myr)
    pub max_erosion_rate_m_per_myr: f32,
    /// Maximum land submergence rate per Myr (m/Myr) 
    pub max_land_submerge_rate_m_per_myr: f32,
    /// Elevation cap minimum (m)
    pub elevation_cap_min_m: f32,
    /// Elevation cap maximum (m)
    pub elevation_cap_max_m: f32,
    /// Small epsilon for numerical comparisons
    pub epsilon: f32,
}

impl Default for PhysConsts {
    fn default() -> Self {
        Self {
            rho_w_kg_per_m3: 1030.0,
            rho_c_kg_per_m3: 2850.0,
            rho_m_kg_per_m3: 3300.0,
            rho_air_kg_per_m3: 1.2,
            g_m_per_s2: 9.81,
            th_ref_continental_m: 35_000.0,
            r_earth_m: 6_371_000.0,
            max_erosion_rate_m_per_myr: 0.5,
            max_land_submerge_rate_m_per_myr: 0.5,
            elevation_cap_min_m: -11_000.0,
            elevation_cap_max_m: 9_000.0,
            epsilon: 1e-6,
        }
    }
}

/// Returns the engine version string from Cargo metadata.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn version_is_semver_like() {
        assert!(version().split('.').count() >= 3);
    }
}
