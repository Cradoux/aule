# PR for {TASK_ID}: {TASK_TITLE}

## Summary
- Task card: {link or ID}
- Scope: {one-liner}

## Changes
- Files touched (only those allowed by the card):
  - …

## Acceptance Criteria Matrix
| Criterion | Evidence |
|---|---|
| C1: {quote from card} | {tests/screens, code refs} |
| C2: {…} | {…} |

## Validation
- Unit tests: `cargo test -p engine` → {pass/fail}
- Lints: `cargo clippy -- -D warnings` → {pass}
- Format: `cargo fmt -- --check` → {pass}
- Viewer smoke test: `cargo run -p viewer` (60 FPS idle) → {ok}

## Benchmarks (include numbers)
- Machine/GPU: {e.g., 13700K + RTX 4080}
- Grid: F=256 (cells≈), Step time: {ms}
- Grid: F=512, Step time: {ms}
- Notes: {any perf regressions/mitigations}

## Screenshots / Plots (if applicable)
- Plate overlay / boundaries / age–depth plot

## Docs
- Updated: README, docs pages, public API rustdoc

## Risk/Assumptions
- {any assumptions taken to resolve ambiguity}

## Checklist
- [ ] Matches card ID & title
- [ ] Only allowed files modified
- [ ] New/changed APIs documented
- [ ] Tests updated
- [ ] Benchmarks run & recorded
- [ ] CI green
- [ ] Determinism maintained