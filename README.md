# Aulë

[![Docs](https://img.shields.io/badge/docs-mdBook-blue)](https://cradoux.github.io/aule/)

Quick start

- Install Rust (stable): see `https://rustup.rs`
- Build workspace:
  - Format check: `cargo fmt -- --check`
  - Lint: `cargo clippy -- -D warnings`
  - Tests: `cargo test --all`
- Run the viewer:
  - `cargo run -p viewer`

Docs

- Browse the docs site: https://cradoux.github.io/aule/
- Build locally:
  - `cargo install mdbook --locked`
  - `mdbook serve docs -p 3000`

Notes

- The viewer opens a window and targets ~60 FPS idle via vsync.
- Engine and export crates are stubs for now; Python crate is a placeholder with optional bindings not enabled by default.

See CONTRIBUTING.md for coding standards and PR flow.

