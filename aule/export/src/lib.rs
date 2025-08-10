//! Export crate stub.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

/// Placeholder function.
pub fn placeholder() -> bool {
    true
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert!(placeholder());
    }
}
