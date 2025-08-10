//! Minimal Python placeholder crate without bindings (to be added in T-201).

/// Returns the package version from Cargo metadata.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn version_present() {
        assert!(version().split('.').count() >= 3);
    }
}
