# Aulë Backlog — Simulation & Architecture Tickets (for Cursor)

_You are **Cursor** (junior dev). I’m the reviewer. Work ticket-by-ticket. For every PR: add tests, update docs, and paste a before/after diagnostic snapshot (first 3 lines of logs + one PNG). Use clear commit messages. Units must be explicit in code and PR text._

---

## Epic 1 — Rate Laws & Numerical Stability

### RL-1 — Unified CFL limiter based on physical widths (replace per-module ad‑hoc substeps)
**Context**: Logs show repeated CFL violations and automatic jumps to 32× substeps (e.g., transforms `v_t*dt/w_half ≈ 1.9`). Substepping hides an overly aggressive update.

**Spec**:
- For any process that moves material a distance `Δ = u * dt` across a characteristic width `W`, compute `C = |Δ| / W` and scale the update by `min(1, CFL_MAX / C)` with `CFL_MAX = 0.3` (tunable).
- Keep substeps only as a last-resort fallback (debug flag).

**Deliverables**:
- `cfl::limit(raw_displacement_m, width_m, dt_myr, cfl_max) -> scaled_displacement_m` utility.
- All modules (transforms, subduction/orogeny, rifting, accretion) call this once per update.
- Log per-process mean/max CFL each step.

**Acceptance**:
- No automatic substep escalations during a 10 Myr run at default dt.
- `max(C)` ≤ `CFL_MAX + 5%` in logs.

---

### RL-2 — Replace hard caps with smooth saturators (on **rates**, not accumulated Δ)
**Context**: We clamp rifting to 150 m/Myr and orogeny to 300 m/Myr almost every step.

**Spec**:
- Implement `soft_cap(x, cap) = cap * tanh(x / cap)`.
- Apply to per-process **rate** before integrating; remove or greatly raise the old hard caps (keep a debug assertion).

**Acceptance**:
- Cap-hit percentage per process < 1% over 10 Myr (reported by new diagnostics).

---

### RL-3 — Derive rates from kinematics with characteristic widths
**Context**: Rates are too large because widths are implicit or too small (e.g., `w_half = 10 km`).

**Spec**:
- **Transforms**: shear strain rate `γ̇ ≈ v_t / W_t`. Use base `W_t = 40–60 km`, or `W_t = clamp(0.5 * litho_thickness, 20 km, 80 km)`.
- **Compression/orogeny**: vertical thickening `ḣ_c = α * v_n * sin(θ)^β / W_c` with `W_c = 150–300 km`, `α ∈ [0.3, 1]`, `β = 1.5`.
- **Rifting**: thinning tied to extension `ḣ_c = -(v_rift / W_r) * h_c` with `W_r = 60–120 km` **or** McKenzie-style β(t) with relaxation `τ = 5–15 Myr` (`dβ/dt = (β* - β)/τ`).
- **Oceanization**: after RL‑3 is in, lower `k_c_oceanize` toward `0.01–0.02 / Myr` and let kinematics dominate.

**Acceptance**:
- Typical rift thinning rates 10–150 m/Myr; orogeny thickening 10–100 m/Myr; transform CFL below 0.3 without substeps.
- Land fraction trend increases (not required to hit a target yet).

---

### RL-4 — Remove per-module substep spam; single scheduler scales updates
**Spec**:
- A small scheduler queries each process for `(raw_displacement, width)` and applies RL‑1 scaling centrally in a fixed order.
- Keep a debug flag to show what scaling was applied per process.

**Acceptance**:
- Logs show one consolidated CFL summary per step; no repeated “increasing substeps …” lines.

---

## Epic 2 — Boundary Logic & Process Exclusivity

### BL-1 — Exclusivity masks (ridge → subduction → transform)
**Spec**:
- Compute per-edge/per-cell `bitflags` of eligible processes.
- Resolve once per step using priority: Ridge (spreading) > Subduction > Transform. Once a flag is consumed, lower-priority processes cannot modify the same cell/edge this step.

**Acceptance**:
- Zero occurrences of stacked edits to the same cell in a single step (assert & log counter).

---

### BL-2 — Subduction obliquity weighting
**Spec**:
- Multiply effective convergence by `sin(θ)^β` with `β = 1.5`.
- Apply consistently across slab consumption, arc/back‑arc budgets, and crustal thickening.
- Keep τ_conv gate (min ≥ 0.020 m/yr; viewer default 0.025) and make it a parameter.

**Acceptance**:
- Oblique segments (θ small) rarely pass τ_conv + obliquity product; arc/back-arc per‑Myr edits drop on oblique margins without relying on caps.

---

### BL-3 — Boundary-prox buoyancy scaling (short-term tectonic dominance bands)
**Spec**:
- Compute distance-to-boundary `d(x)` (cells to nearest active ridge/trench/transform).
- Weight isostatic response with a smooth kernel: `w_b(d) = 0.5 * [1 + tanh((d - d0)/s)]`.
- Defaults: `d0 = 50 km`, `s = 20 km`. Multiply vertical buoyancy update by `w_b(d)`.

**Acceptance**:
- Trust-region hits in buoyancy (if any left) occur **away** from active boundaries; vertical noise bands at boundaries reduce visibly in PNGs.

---

## Epic 3 — Isostasy/Buoyancy Solver

### IB-1 — Implicit or line‑searched isostatic update (remove 200 m trust region)
**Spec**:
- Option A: solve vertical equilibrium implicitly for `Z^{k+1}` given loads/densities.
- Option B: keep explicit update but wrap with backtracking line search (Armijo) on `ΔZ` to ensure energy decrease; remove the blunt `±200 m` per-step clamp.

**Acceptance**:
- Logs: `dZ_step max` is no longer pinned at a fixed cap; no oscillatory overshoots; “amp” stabilizes under ~1 km for typical loads.

---

### IB-2 — Density & unit audit
**Spec**:
- Centralize crust/mantle/water/air densities; ensure consistent units and conversions; document them in `SimParams`.

**Acceptance**:
- One source of truth for densities; a unit test that a neutral column yields `ΔZ=0`.

---

## Epic 4 — Diagnostics & Logging

### DL-1 — Boundary composition & cap metrics
**Spec**:
- Per step: report % boundary length by class; counts of oceanized cores; and % of cells/edges hitting any cap (pre‑ and post‑RL‑2 numbers).

**Acceptance**:
- New section in logs with a one-line summary every Myr; CSV export for viewer plot.

---

### DL-2 — Mass/volume budgets
**Spec**:
- Track net changes in crustal thickness, subducted volume, accreted volume, and vertical integration (isostatic compensation). Verify global conservation within tolerance.

**Acceptance**:
- Budget closes within a small % tolerance (set and logged); failing the tolerance raises a warning.

---

### DL-3 — `tracing` spans & throttled summaries
**Spec**:
- Use `tracing` to group logs per process; throttle repetitive lines; emit per‑Myr summaries with mean/stddev of key metrics.

**Acceptance**:
- Cleaner logs; repeated transform diagnostics collapse into a single summary per step.

---

## Epic 5 — Parameters, Units, and Scheduler Architecture

### AP-1 — Typed `SimParams` (single source of truth)
**Spec**:
- Move all magic numbers (e.g., 0.025, 10_000 m, 150 m/Myr) into `SimParams` with units in field names and doc comments.
- Provide serde (toml) load + profile selection.

**Acceptance**:
- No hard-coded tunables in modules; profiles: `conservative`, `standard`, `aggressive` loadable at startup.

---

### AP-2 — Units-of-measure wrappers
**Spec**:
- Introduce minimal newtypes or adopt a UoM crate for `Meters`, `Kilometers`, `Myr`, `MetersPerYear`.
- Conversions are explicit.

**Acceptance**:
- Failing compile if a `km` value is fed where `m` is expected (via newtypes or `From` conversions).

---

### AP-3 — Central scheduler (single place that orders processes and applies CFL scaling)
**Spec**:
- Scheduler requests `(raw_displacement, width)` from each process, applies RL‑1, then commits updates in fixed order: Ridge → Subduction → Transform → Accretion/Orogeny → Isostasy/Flexure.

**Acceptance**:
- One unified step loop; modules no longer self‑manage substeps.

---

### AP-4 — Numerical consistency
**Spec**:
- Use `f64` for accumulation; keep storage `f32` where required; avoid re-building local tangents—route through one geo utility.

**Acceptance**:
- Single source of truth for local basis; no duplicate tangent calculators.

---

## Epic 6 — Tests

### TST-1 — Isostasy property tests
**Spec**:
- A neutral reference column returns `ΔZ ≈ 0`.
- Doubling crustal thickness increases elevation consistent with densities (tolerance specified).

**Acceptance**:
- Tests pass; tolerance documented.

---

### TST-2 — Invariants & guards
**Spec**:
- Assert no negative crustal thickness; slopes within reasonable bounds; exclusivity prevents stacked edits.

**Acceptance**:
- Failing invariants panic in debug mode and log in release with counts.

---

### TST-3 — Snapshot step on tiny grid (8×4)
**Spec**:
- One step snapshot test stores: boundary counts, land %, cap‑hit %, and CFL stats.

**Acceptance**:
- Deterministic within a tolerance; used in CI to catch regressions.

---

### TST-4 — CFL limiter unit tests
**Spec**:
- Given `(raw, W, dt)`, ensure `|raw_scaled|/W ≤ CFL_MAX` and monotonicity.

**Acceptance**:
- Tests cover normal & edge cases.

---

## Epic 7 — Viewer & UX Diagnostics

### UI-1 — On‑screen diagnostics
**Spec**:
- Display boundary-length percentages, oceanized counts, cap‑hit %, and land % as overlays/legend blocks.

**Acceptance**:
- Toggling debug overlay shows these values updating each step.

---

### UI-2 — Parameter profiles & tooltips
**Spec**:
- Expose profile selector; show unit‑annotated tooltips for major knobs (dt, widths, caps, τ_conv, β, α, CFL_MAX).

**Acceptance**:
- Switching profiles updates `SimParams` live (with confirmation prompt).

---

### UI-3 — Hypsometry & land‑fraction plots vs time
**Spec**:
- Simple time‑series panel plotting land %, mean elevation, and cap‑hit % from DL‑1 CSV.

**Acceptance**:
- 60 Myr run shows trend; exported PNG/CSV via UI.

---

## How to implement (order of work)
1. **BL‑1 (Exclusivity)** and **BL‑2 (Obliquity)** — low risk, immediate stability gains.
2. **RL‑1 (CFL limiter)** + **RL‑4 (Scheduler)** — removes substep spam.
3. **RL‑2 (Smooth caps)** + **RL‑3 (Width‑based laws)** — makes rates honest.
4. **IB‑1 (Isostasy line search/implicit)** + **BL‑3 (Buoyancy scaling)** — removes trust region.
5. **Diagnostics (DL‑1..3)** and **Tests (TST‑1..4)** — lock in behavior.
6. **AP‑1..4** and **UI‑1..3** as follow‑ups to tidy architecture & UX.

---

## Definition of Done (for each ticket)
- Code + tests + docs + param updates.
- Logs include before/after metrics.
- No new warnings; unit consistency checks pass.
- Reviewer (me) signs off with a short benchmark (10 Myr, default dt) PNG and log snippet.

---

## Notes & Starting Values (tunable)
- `CFL_MAX = 0.30`.
- `W_t = 50 km` (min 20, max 80), `W_c = 200 km`, `W_r = 90 km` to start.
- `β_obliquity = 1.5`, `α_orogeny = 0.5`.
- Buoyancy boundary scaling: `d0 = 50 km`, `s = 20 km`.
- Revisit `k_c_oceanize` after RL‑3; target `0.01–0.02 / Myr`.

> Cursor: start with **BL‑1**, **BL‑2**, then **RL‑1**. Ping me with a PR after each, including a 2D map PNG + the single‑step CFL summary line from the logs.

