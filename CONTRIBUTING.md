# Contributing to Aulë

## Local dev loop
1) Format:        cargo fmt --all
2) Lints (strict): cargo clippy --all-targets -- -D warnings
3) Tests:         cargo test --all
4) Run viewer:    cargo run -p viewer

## CI parity
The CI runs the same three commands as above (fmt/clippy/test) on Linux/Windows/macOS. Keep them green.

## Branch & PR flow
- Branch name: T-XXX-short-title (e.g., T-010-grid)
- One PR per task card. Scope = card’s Deliverables only.
- Use .github/pull_request_template.md and fill the Acceptance Criteria Matrix.
- Include OS/GPU + any perf numbers if relevant.

## Coding standards
- No `unsafe` in MVP.
- No `unwrap`/`expect`/`dbg!` in non-test code.
- Public APIs documented (especially in `engine`).
- Deterministic runs for a fixed seed.

## Dependencies
- Keep to the crates approved in the task cards. Ask before adding anything else.

## Definition of Done
- Acceptance criteria met, tests/docs updated, CI green, no TODOs left in code.
