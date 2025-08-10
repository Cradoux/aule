//! Aulë engine crate stub.
//! Minimal, no GPU yet.
#![deny(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

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
