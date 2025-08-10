//! Aulë engine crate stub.
//! Minimal; includes GPU helpers for tests.
#![deny(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

/// Field and tiling views.
pub mod fields;
<<<<<<< HEAD
=======
/// Geodesic grid module.
pub mod grid;
/// Plates seeding and velocities.
pub mod plates;
>>>>>>> cb7bf80 (T-030: plate seeds (FPS), Voronoi labeling, Euler poles, per-cell velocities; viewer logs |V| stats at F=64)
/// Minimal GPU helper for tests (T-020).
pub mod gpu;
/// Geodesic grid module.
pub mod grid;

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
