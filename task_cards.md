# Aulë – Cursor Task Cards (v0.1)

> **Source of truth:** Follow these cards exactly. If anything is ambiguous, ask the Architect (lead model) before deviating. Do not introduce new dependencies without an explicit card.

---

## Task template (use for every PR)

**ID & Title**
**Objective**: *Single sentence outcome.*
**Deliverables**: *Files created/modified, outputs, docs.*
**Implementation outline**: *Bullet steps only.*
**Acceptance criteria**: *Binary checks + perf budget.*
**Dependencies**: *IDs.*
**Notes**: *Assumptions, guardrails.*

**Definition of Done (DoD)**

* All acceptance criteria met.
* Unit tests added/updated.
* Benchmarks (if perf‑related) updated and pass on CI.
* Docs: README section or docs page updated.
* No dead code, no TODOs left in code.

---

## Sprint 0 – Scaffolding (Week 1)

### T‑000 — Repo bootstrap & CI

**Objective:** Initialize Rust + wgpu project with viewer and engine crates; set up CI.
**Deliverables:** `/aule/engine`, `/aule/viewer`, `/aule/python`, `/docs` skeleton; GitHub Actions (build/test on Win/Linux/macOS); `LICENSE` (Apache‑2.0), `README.md` with run instructions; `rust-toolchain.toml`.
**Implementation outline:**

* Create Cargo workspace: `engine` (lib), `viewer` (bin), `export` (lib), `python` (pyo3 placeholder).
* Add `wgpu`, `winit`, `egui`, `thiserror`, `bytemuck`, `glam`, `rand`.
* CI: cache Rust; run clippy + fmt + tests.
  **Acceptance criteria:** CI green on all OSes; `cargo run -p viewer` opens empty window @ 60 FPS.
  **Dependencies:** None.
  **Notes:** No GPU kernels yet; just device setup and hot‑reload loop.

### T‑005 — Coding standards & lints

**Objective:** Enforce consistent style and safe patterns.
**Deliverables:** `rustfmt.toml`, `clippy.toml`, `CONTRIBUTING.md`, PR template.
**Implementation outline:** configure strict clippy (deny warnings), forbid `unsafe` in `engine` for now.
**Acceptance criteria:** `cargo clippy -- -D warnings` passes.
**Dependencies:** T‑000.

### T‑008 — Docs site skeleton

**Objective:** Minimal docs with mkdocs or mdBook.
**Deliverables:** `/docs` with index, build on CI to GitHub Pages (if available).
**Acceptance criteria:** Site builds; links to README and architecture diagram.
**Dependencies:** T‑000.

---

## Sprint 1 – Grid & Data (Week 1)

### T‑010 — Geodesic grid & adjacency

**Objective:** Implement icosphere grid (`F` frequency) with cell‑centered dual mesh and neighbor lists.
**Deliverables:** `engine/src/grid.rs`; binary cache writer/reader for grid.
**Implementation outline:**

* Build icosahedron → subdivide → project to sphere.
* Compute dual cell centers (Voronoi of vertices) and store 1‑ring/2‑ring adjacency.
* Expose `Grid { cells, neighbors, area, latlon }`.
  **Acceptance criteria:** Unit tests: uniform area within 5%; neighbor counts stable; serialization round‑trips.
  **Dependencies:** T‑000.

### T‑015 — Tile partitioning

**Objective:** Partition grid into tiles for GPU dispatch.
**Deliverables:** `engine/src/grid.rs` (tiling API), `engine/src/fields.rs` (tile views).
**Implementation outline:** Equal‑size tiles with halos; mapping to WGSL dispatch groups.
**Acceptance criteria:** Tiles cover all cells; halos complete; test for boundary wrap correctness.
**Dependencies:** T‑010.

### T‑020 — GPU field buffers (SoA)

**Objective:** Allocate structured device buffers for core fields.
**Deliverables:** `engine/src/fields.rs`, `shaders/buffers.wgsl`.
**Implementation outline:** SoA buffers (f32 for scalars, u16 for IDs), staging for readback; zero‑copy slices.
**Acceptance criteria:** Create/resize at runtime; round‑trip write/read fuzz test.
**Dependencies:** T‑015.

---

## Sprint 2 – Kinematics & Boundaries (Week 2)

### T‑030 — Plate seeds & Euler poles

**Objective:** Seed N plates; compute per‑cell velocities from Euler poles.
**Deliverables:** `engine/src/plates.rs`.
**Implementation outline:** Farthest‑point sampling for seeds → Voronoi labels; random seeded Euler poles; compute surface velocity `V` per cell.
**Acceptance criteria:** Deterministic output given `seed`; visual sanity (viewer overlay) shows arrows and labels.
**Dependencies:** T‑020.

### T‑040 — Boundary detection

**Objective:** Classify divergent/transform/convergent boundaries.
**Deliverables:** `engine/src/boundaries.rs`.
**Implementation outline:** From neighbor plate diffs + relative velocity decomposition; produce boundary type field `B`.
**Acceptance criteria:** Unit test on synthetic 3‑plate setup; correct counts/types.
**Dependencies:** T‑030.

### T‑120a — Viewer overlays (plates)

**Objective:** Visualize plate IDs, velocity vectors, and boundary types.
**Deliverables:** `viewer`: layers & legend; key toggles.
**Acceptance criteria:** Hotkeys: `1` plates, `2` velocities, `3` boundaries.
**Dependencies:** T‑040.

---

## Sprint 3 – Oceans: Ridges, Age, Bathymetry (Week 3)

### T‑050 — Ridge kernel

**Objective:** Create zero‑age oceanic crust at divergent boundaries.
**Deliverables:** `engine/src/boundaries.rs`, `shaders/ridge.wgsl`.
**Implementation outline:** Identify ridge cells; reset `age_ocean=0`; initialize oceanic thickness proxy; smear to 1‑ring.
**Acceptance criteria:** Age histogram shows births at ridges; land fraction remains stable.
**Dependencies:** T‑040, T‑020.

### T‑060 — Oceanic age & age→depth

**Objective:** Advect `age_ocean`; compute bathymetry from age.
**Deliverables:** `engine/src/bathymetry.rs`, `shaders/age_advect.wgsl`.
**Implementation outline:** Increment age by Δt off‑ridge; depth = piecewise function of age; apply water‑loaded correction.
**Acceptance criteria:** Age–depth curve within tolerance vs reference; visual abyssal deepening.
**Dependencies:** T‑050.

### T‑061 — Age–depth validation plot

**Objective:** Diagnostic plot in viewer.
**Deliverables:** viewer chart module.
**Acceptance criteria:** Overlay sim vs reference; RMS misfit threshold.
**Dependencies:** T‑060.

---

## Sprint 4 – Convergence: Subduction & Arcs (Week 4–5)

### T‑070 — Subduction kernel

**Objective:** Consume oceanic crust; create trench + arc uplift.
**Deliverables:** `engine/src/boundaries.rs`, `shaders/subduction.wgsl`.
**Implementation outline:** Identify slab‑hinge at convergent cells; deepen trench; uplift fore‑arc/arc/back‑arc bands; remove oldest ocean floor.
**Acceptance criteria:** Trench and arc masks export; land/ocean % reasonable; no negative ages.
**Dependencies:** T‑040, T‑060.

### T‑071 — Rollback & back‑arc extension

**Objective:** Add optional trench rollback and back‑arc thinning.
**Deliverables:** params + kernel branch; config toggles.
**Acceptance criteria:** Back‑arc basins appear under rollback>0; bounded by stability tests.
**Dependencies:** T‑070.

### T‑072 — Transforms (continental & oceanic)

**Objective:** Strike‑slip offsets and pull‑apart basins.
**Deliverables:** kernel rules; masks.
**Acceptance criteria:** Visible offsets; basin statistics within limits.
**Dependencies:** T‑040.

---

## Sprint 5 – Flexure & Isostasy (Week 5–6)

### T‑080 — Flexural solver (variable Te)

**Objective:** Thin‑plate elastic flexure with multigrid on tiles.
**Deliverables:** `engine/src/flexure.rs`, `shaders/flexure_*.wgsl`.
**Implementation outline:** Assemble loads (topo+sediment+water); 2–3 V‑cycles per step; mixed BC via tiling halos.
**Acceptance criteria:** Line‑load analytic test within 5%; converges in < 10 iterations typical.
**Dependencies:** T‑020, T‑100 (sediment load later).

### T‑090 — Airy coupling

**Objective:** Local buoyancy snap without double counting against flexure.
**Deliverables:** `engine/src/isostasy.rs`.
**Acceptance criteria:** No drift in mean sea level; hills/mountains respond sensibly to erosion changes.
**Dependencies:** T‑080.

---

## Sprint 6 – Surface Processes (Week 6–7)

### T‑100 — Hydrology & erosion (stream‑power)

**Objective:** Flow routing + incision + diffusion.
**Deliverables:** `engine/src/hydro.rs`, `shaders/hydro_*.wgsl`.
**Implementation outline:** MFD routing on sphere; stream‑power `E = K A^m S^n`; hillslope diffusion; coastal flooding & retreat.
**Acceptance criteria:** Mass conservation test; stable hypsometry; rivers reach sea.
**Dependencies:** T‑020.

### T‑101 — Sediment routing & deposition

**Objective:** Exner mass balance; shelves/basins fill.
**Deliverables:** same modules as T‑100.
**Acceptance criteria:** Sediment thickness non‑negative; foreland basins form under loads.
**Dependencies:** T‑100, T‑080.

### T‑102 — Surface process benchmarks

**Objective:** Perf tests and stability sweeps.
**Deliverables:** `benchmarks/` + CI job.
**Acceptance criteria:** F=256 step < 100 ms on target GPU; no blow‑ups in 1000‑step sweeps.
**Dependencies:** T‑100.

---

## Sprint 7 – Climate (Week 7)

### T‑110 — Latitudinal winds & base precip

**Objective:** Zonal winds (Hadley/Ferrel/Polar) + base precip field.
**Deliverables:** `engine/src/climate.rs`.
**Acceptance criteria:** Precip peaks near ITCZ; polar deserts; viewer map.
**Dependencies:** T‑020.

### T‑111 — Orographic precipitation (linear)

**Objective:** Orographic enhancement using linear model; optional monsoon factor.
**Deliverables:** shader or CPU helper; config toggles.
**Acceptance criteria:** Rain shadows on leeward sides; precip histograms look plausible.
**Dependencies:** T‑110.

---

## Sprint 8 – Viewer, Exports, Polish (Week 8)

### T‑120 — Viewer controls & diagnostics

**Objective:** Sliders & plots.
**Deliverables:** toggles/knobs (plates, ridge strength, Te, rift threshold, erosion K, monsoon, plume rate, SLR, seed); charts (hypsometry, age–depth, land %, precip histogram).
**Acceptance criteria:** Real‑time updates; screenshot/export works.
**Dependencies:** previous viewer subtasks.

### T‑130 — Exporters

**Objective:** PNG16/TIFF32/NetCDF or Zarr + masks + plate vectors (GeoJSON/GPML).
**Deliverables:** `/export` crate; examples.
**Acceptance criteria:** Files readable in Blender/QGIS/Gaea; metadata included; legends.
**Dependencies:** T‑020.

### T‑131 — Color atlases & legends

**Objective:** Hypsometric tint, shaded relief preview, biome preview.
**Deliverables:** `/assets/` ramps; preview shaders.
**Acceptance criteria:** Consistent colors across exports and viewer.
**Dependencies:** T‑130.

---

## Optional – Integrations & CLI

### T‑140 — pyGPlates bridge (kinematic mode)

**Objective:** Import plate polygons + rotation model; reconstruct to grid.
**Deliverables:** `python/` bindings; `scripts/gplates_import.py`.
**Acceptance criteria:** Known GPML sample reconstructs; plate IDs/velocities match.
**Dependencies:** T‑030, T‑040.

### T‑200 — CLI runner

**Objective:** Headless runs from config.
**Deliverables:** `run.py`; `configs/` samples (earthlike, archipelago, supercontinent).
**Acceptance criteria:** Reproducible outputs with seed; logs with timings.
**Dependencies:** Core engine ready.

### T‑201 — Python bindings

**Objective:** pyo3 wrapper for core APIs.
**Deliverables:** `python/lib.rs`, wheel build on CI.
**Acceptance criteria:** `pip install` local wheel; minimal example.
**Dependencies:** T‑200.

---

## Cross‑cutting

### T‑300 — Benchmarks

**Objective:** Micro + end‑to‑end benchmarks.
**Acceptance criteria:** Baselines checked in; regressions fail CI.

### T‑310 — Deterministic RNG

**Objective:** Seeded RNG for stochastic events.
**Acceptance criteria:** Same inputs → same outputs.

### T‑320 — Logging & telemetry

**Objective:** Structured logs; perf counters per kernel.
**Acceptance criteria:** Viewer HUD shows timings; JSON logs written.

---

## Guardrails (don’t deviate)

* Stay within Rust + wgpu + minimal deps.
* No external heavy geodynamics libs in MVP.
* Do not change file layout without card.
* Ask before adding crates with transitive GPU/FFT dependencies; prioritize in‑house WGSL or simple Rust.
* Performance budgets are hard limits in acceptance criteria.

## Non‑goals (MVP)

* Full mantle convection, visco‑plastic rheology, or 3D thermal evolution.
* Perfect geologic time scales; we operate in stylized Myr steps.
* GIS projections beyond simple exports.

---

## Milestone exit criteria

**MVP demo**: interactive viewer, plate drift + spreading, subduction shaping, flexure response, basic erosion & climate, exports + masks, validation plots; 500 Myr @ F=512 runs headless in < 10 minutes on target GPU.

---

## Attachments for Cursor setup (paste-ready)

### `CURSOR_SYSTEM_PROMPT.md`

**Paste this entire block into Cursor → Workspace Settings → Instructions. Also commit it at repo root as `CURSOR_SYSTEM_PROMPT.md`.**

```
Title: Aulë Implementer Agent — System Prompt

ROLE
You are the Implementer agent for the Aulë project (Procedural Tectonics & World Builder).

SOURCE OF TRUTH
1) "Tectonics GPU Simulator – Technical Blueprint v0.1" (architecture, science, kernels).
2) "Aulë – Cursor Task Cards v0.1" (this repo’s authoritative backlog).
If a conflict exists, ask the Architect. Do not invent requirements.

OPERATING PRINCIPLES
- Only work on the currently assigned task card.
- Touch only files listed under Deliverables for that card.
- Preserve determinism: fixed seed ⇒ identical outputs.
- Keep performance budgets; if not met, stop and ask.
- No new dependencies/large crates/FFTs unless the card says so.
- Prefer Rust + wgpu (WGSL) with SoA buffers; avoid unsafe code.
- All public APIs documented; tests and benchmarks updated.

WORKFLOW
1) PLAN: Restate the card in your own words; list files to change; list acceptance criteria as checkboxes.
2) IMPLEMENT: Write minimal code to satisfy criteria; keep commits small and readable.
3) VALIDATE: Run unit tests, viewer smoke test, and benchmarks specified by the card. Collect timings.
4) REPORT: Prepare a PR that maps each change to each acceptance criterion. Include step timings and screenshots when relevant.

ASK POLICY
- If ambiguous, ask up to 3 crisp questions in one message. If unanswered, default to the simplest choice that satisfies the card and clearly state the assumption in the PR.

DETERMINISM & LOGGING
- All stochastic paths must be seeded via the project RNG.
- Emit structured perf logs per kernel.

GUARDRAILS
- Do not refactor outside the card’s scope.
- Do not change file layout or add deps without a card.
- Keep CI green (fmt, clippy with -D warnings, tests, benches).

DEFINITION OF DONE (DoD)
- Acceptance criteria met, tests added/updated, benchmarks updated, docs updated, CI green. No TODO/FIXME left in code.

PR FORMAT
- Use the provided PR template. Include an Acceptance Criteria Matrix.

CHECKLIST BEFORE OPENING A PR
- [ ] Matches task ID & title
- [ ] Only allowed files modified
- [ ] New/changed APIs documented
- [ ] Unit tests cover new logic
- [ ] Benchmarks run & recorded
- [ ] Viewer runs @ ~60 FPS idle
- [ ] Deterministic outputs preserved
```

---

### `.github/pull_request_template.md`

**Create this file to standardize PRs.**

```
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
```

---

### `CHECKS.md`

**Optional quick list Cursor can run through before each PR.**

```
- Task matches current card; no scope creep
- No new dependencies added
- Deterministic seed path intact
- Perf budgets met or flagged
- Exports unchanged unless card says so
- Logs include kernel timings
- Viewer hotkeys (1/2/3) work after change
```

---

### Suggested repo files updated by T‑000

Add stubs so Cursor has something concrete to fill:

```
/README.md                 # quick start, build/run
/CURSOR_SYSTEM_PROMPT.md   # paste into Cursor and commit here
/.github/pull_request_template.md
/CHECKS.md
```

> When you assign the next card to Cursor, paste: “Work on {TASK\_ID}. Use CURSOR\_SYSTEM\_PROMPT.md and follow the card strictly. When done, open a PR with the template.”
