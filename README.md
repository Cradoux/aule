# Aulë

Quick start

- Install Rust (stable): see `https://rustup.rs`
- Build workspace:
  - Format check: `cargo fmt -- --check`
  - Lint: `cargo clippy -- -D warnings`
  - Tests: `cargo test --all`
- Run the viewer:
  - `cargo run -p viewer`

Notes

- The viewer opens a window and targets ~60 FPS idle via vsync.
- Engine and export crates are stubs for now; Python crate is a placeholder with optional bindings not enabled by default.

