//! Aulë engine crate stub.
//! Minimal, no GPU yet.

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
