# Aulë – Cursor Task Cards (v0.2)

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

## GitHub workflow (for every task)

**Branching policy**

* Start from up-to-date `main`.
* Create a new branch named with the pattern: `{type}/T-XXX-kebab-title` where `{type}` ∈ {`feat`, `fix`, `chore`, `docs`}.
* Keep scope to the card’s Deliverables only. If a task depends on another in-flight task, open as a **Draft PR** and mark the dependency in the PR description.

**Commands (copy/paste)**

*PowerShell*

```powershell
# ensure main is current
git checkout main; git pull
# create task branch (example for T-015)
git checkout -b feat/T-015-tile-partitioning
# work…
git add -A; git commit -m "feat(T-015): implement tile partitioning"
# publish branch and open PR
git push -u origin HEAD
```

*Bash*

```bash
# ensure main is current
git checkout main && git pull
# create task branch (example for T-015)
git checkout -b feat/T-015-tile-partitioning
# work…
git add -A && git commit -m "feat(T-015): implement tile partitioning"
# publish branch and open PR
git push -u origin HEAD
```

**PR rules**

* One PR per task. Use the template. Title: `T-XXX — {Title}`.
* CI must be green (fmt, clippy, tests). If blocked by another PR, keep as **Draft**.
* Prefer **Squash & merge**; delete the branch after merge.

**Keeping in sync**

*PowerShell*

```powershell
git fetch origin
# rebase your branch on latest main (preferred)
git rebase origin/main
# or merge if necessary (avoid unless conflicts are complex)
# git merge origin/main
```

*Bash*

```bash
git fetch origin
# rebase your branch on latest main (preferred)
git rebase origin/main
# or merge if necessary (avoid unless conflicts are complex)
# git merge origin/main
```

\*\*Branch types by card\*\*

\- T-000: \`feat/T-000-bootstrap\`

\- T-005: \`chore/T-005-lints\`

\- T-008: \`docs/T-008-docs-site\`

\- T-010: \`feat/T-010-grid\`

\- T-015: \`feat/T-015-tiling\`

\- T-020: \`feat/T-020-fields\` (\*\*Draft\*\* until T-015 merges)

\- T-030: \`feat/T-030-plates\`

\- T-040: \`feat/T-040-boundaries\`

\- T-051: \`feat/T-051-continents\`

\- T-060: \`feat/T-060-age-bathy\`

\- T-061: \`feat/T-061-age-depth-plot\`

\- T-070: \`feat/T-070-subduction\`

\- T-071: \`feat/T-071-rollback\`

\- T-072: \`feat/T-072-transforms\`

\- T-073: \`feat/T-073-edge-kinematics\`

\- T-080a: \`feat/T-080a-flexure-ref\`

\- T-080b: \`feat/T-080b-flexure-gpu\`

\- T-082: \`feat/T-082-flexure-coupling\`

\- T-090: \`feat/T-090-isostasy\`

\- T-100: \`feat/T-100-hydro-erosion\`

\- T-101: \`feat/T-101-sediment\`

\- T-102: \`chore/T-102-benchmarks\`

\- T-110: \`feat/T-110-winds-precip\`

\- T-111: \`feat/T-111-orographic\`

\- T-120: \`feat/T-120-viewer-controls\`

\- T-130: \`feat/T-130-exporters\`

\- T-131: \`feat/T-131-color-atlases\`

\- T-140: \`feat/T-140-pygplates\`

\- T-200: \`feat/T-200-cli\`

\- T-201: \`feat/T-201-python-bindings\`

\- T-300: \`chore/T-300-benchmarks\`

\- T-310: \`chore/T-310-deterministic-rng\`

\- T-320: \`chore/T-320-logging\`

**Conventional commits (recommended)**

* `feat(T-015): add tiler with halo=2`
* `fix(T-010): correct UP loop bounds`
* `docs(T-008): add site-url=/aule/ and README badge`
* `chore(T-005): tighten clippy and add CONTRIBUTING`

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

\### T-073 — Persist per-edge kinematics

\*\*Objective:\*\* Store \`(n, t, n̂, t̂)\` once during boundary classification and reuse in subduction/transforms.

\*\*Deliverables:\*\* Extend \`Boundaries\` with \`edge\_kin: Vec\<EdgeKin>\`; refactor \`subduction.rs\` and \`transforms.rs\` to consume; tests comparing stored vs recomputed \`(n,t)\`; no behavior change.

\*\*Implementation outline:\*\* Canonicalize edges \`u\<v\`; compute midpoint frame and project velocities via shared helper; sort by \`(u,v)\`.

\*\*Acceptance criteria:\*\* Stored vs ad-hoc \`(n,t)\` agree within 1e-9; subduction/transforms masks & stats unchanged; CI green.

\*\*Dependencies:\*\* T-040.

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

\### T-051 — Continental mask (synthetic)

\*\*Objective:\*\* Deterministic large-scale cratons to create land; overlay + amplitude targeting.

\*\*Deliverables:\*\* \`engine/src/continent.rs\` (caps + smooth union + template, amplitude solve & apply), unit tests; viewer HUD & overlay (toggle \*\*C\*\*) with seed/n/radius/falloff, auto-amplitude to target land %, manual amplitude, cached meshes; recompute order before sea level.

\*\*Implementation outline:\*\*

\- Place \`n\` spherical caps from seed; plateau inside mean radius; Gaussian falloff; smooth-union into \[0,1] template.

\- Solve amplitude by bisection to hit target area-weighted land fraction, else use manual amplitude.

\- Apply uplift: \`depth += -(amp \* template)\`; build land mask & coastline mesh; cache.

\*\*Acceptance criteria:\*\* Deterministic template; auto target within ±1 pp land %; viewer overlay renders; recompute finishes in a few ms at F=64.

\*\*Dependencies:\*\* T-060 (depth), T-070/T-072 (order in recompute).

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

### ### T-080a — Flexure reference (CPU, 1D)


\*\*Objective:\*\* Analytic line-load reference and CG solver for 1D plate on Winkler foundation; viewer plot (\*\*F\*\*).

\*\*Deliverables:\*\* \`engine/src/flexure.rs\` (D\_from\_Te, alpha, w\_line\_analytic, 1D CG solver), tests (RMS ≤5% vs analytic), viewer plot with iterations/residual/RMS.

\*\*Dependencies:\*\* T-020.

\### T-080b — Flexure WGSL multigrid (tiles)

\*\*Objective:\*\* WGSL kernels and multigrid V-cycle on tiled atlas; manufactured solution residual reduction.

\*\*Deliverables:\*\* \`engine/src/flexure\_gpu.rs\`, \`shaders/flexure.wgsl\` (apply\_A, residual, Jacobi, restrict, prolong), tests: CPU vs GPU operator (1e-5), ignored GPU residual reduction (<0.5).

\*\*Dependencies:\*\* T-020, T-080a.

\### T-082 — Flexure load assembly & viewer coupling

\*\*Objective:\*\* Assemble loads from depth (water/rock), run 1–2 GPU V-cycles, apply deflection in viewer; HUD controls.

\*\*Deliverables:\*\* \`engine/src/flexure\_loads.rs\` (assemble), viewer HUD (G) with params (E,ν,Te,k,ω,ν1,ν2,levels), residual log; optional w overlay.

\*\*Acceptance criteria:\*\* Residual reduction printed; depth changes visibly near loads; toggle on/off works; CI green.

\*\*Dependencies:\*\* T-080b, T-051, T-090.

### ### T-090 — Global sea-level constraint (isostasy MVP)


\*\*Objective:\*\* Maintain constant reference ocean volume by solving a uniform sea-level offset; viewer toggle (\*\*L\*\*) with target fraction.

\*\*Deliverables:\*\* \`engine/src/isostasy.rs\` (ocean volume, bisection solve, apply offset), tests (monotonicity, solver accuracy, idempotence, determinism); viewer HUD (persisted target fraction, extra Δoffset, lockable bathy scale, 0 m coastline).

\*\*Acceptance criteria:\*\* Logs show non-zero offset when hypsometry changes; coastline shift visible with locked scale; CI green.

\*\*Dependencies:\*\* T-060, T-051.

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

**Paste this entire block into Cursor → Workspace Settings → Instructions. Also commit it at repo root as ****************************`CURSOR_SYSTEM_PROMPT.md`****************************.**

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
