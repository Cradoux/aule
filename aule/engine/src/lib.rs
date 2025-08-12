//! Aulë engine crate stub.
//! Minimal; includes GPU helpers for tests.
#![deny(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

/// Steady-state oceanic age and bathymetry mapping.
pub mod age;
/// Plate boundary classification.
pub mod boundaries;
/// Field and tiling views.
pub mod fields;
/// Geometric utilities used by solvers.
pub mod geo;
/// Minimal GPU helper for tests (T-020).
pub mod gpu;
/// Geodesic grid module.
pub mod grid;
/// Plates seeding and velocities.
pub mod plates;
/// Ridge births and fringe assignment (CPU pass).
pub mod ridge;
/// World stepper.
pub mod stepper;
/// Subduction bands and bathymetry adjustments (CPU pass).
pub mod subduction;
/// World state.
pub mod world;

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
