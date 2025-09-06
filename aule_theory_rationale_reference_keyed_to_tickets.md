# Aulé – Theory & Rationale Reference (keyed to tickets)

> This document explains the theory, equations, and suggested implementation patterns behind each ticket in the backlog. It’s written to be unambiguous for Cursor to implement. Ticket keys match the backlog (e.g., **NS‑01**, **SB‑02**). If a key is missing or renamed in the backlog, keep the content and adjust the key in the PR.

---

## Conventions & Units

- **v** = relative plate velocity (m/yr). Decompose into **tangential** (\|t\|) and **normal** (\|n\|) components to the local boundary normal **n̂**.
- **dt** = model step (Myr). Where sub‑stepping is used, **dt_sub = dt / N_sub**.
- **w_half** = half‑width of the process stencil/band (m) (e.g., 10 km for transforms), so full width = 2·w_half.
- **θ** = obliquity = angle between convergence vector and trench normal.
- **th_c**/**th_o** = continental/oceanic lithospheric thickness (m).
- **ρ_c, ρ_o, ρ_m** = densities of crust/oceanic lithosphere/mantle (kg/m³). Typical: ρ_c≈2800, ρ_o≈3000, ρ_m≈3300.
- **CFL** = Courant number.
- Unless noted otherwise: **distances in meters, time in years/Myr, velocities in m/yr**.

---

## NS‑01 — Transform CFL enforcement & dynamic sub‑stepping

**Why**: Stability. We must bound displacement per sub‑step across the active shear band.

**Theory**
- Courant number for shear advection across a half‑width **w_half**:
  
  **CFL_t = (\|t\| · dt) / w_half**.
  
- Require **CFL_t ≤ CFL_max** (use 0.5). Choose **N_sub = ceil(max(1, CFL_t / CFL_max))** and integrate with **dt_sub = dt / N_sub**.

**Pseudocode**
```ts
const CFL_max = 0.5;
const CFL_t = (vt * dtMyr * 1e6) / w_half; // dt in years here
const Nsub = Math.ceil(Math.max(1, CFL_t / CFL_max));
for (let s = 0; s < Nsub; s++) {
  applyTransformShear(dt / Nsub);
}
```

**Notes**
- Log **CFL_t**, **N_sub**, and the top‑k outliers to help tuning.
- Do the same for other band‑based edits (see NS‑02).

---

## NS‑02 — Subduction CFL (banded edits) & consistent sub‑stepping

**Why**: You already observe occasional large CFL for subduction bands.

**Theory**
- Use the same band‑CFL but for the **normal** component driving trench‑perpendicular advection:
  
  **CFL_n = (\|v·n̂\| · dt) / w_half_trench**.

**Implementation**
- Share the utility from NS‑01. Use a potentially larger **w_half_trench** than transforms to reflect wider deformation zones.
- Keep **CFL_max=0.5** unless we demonstrate monotone schemes at higher Courant.

---

## BC‑01 — Exclusivity masks (ridge → subduction → transforms)

**Why**: Prevent stacking of mutually exclusive edits in a single cell per step and enforce geologic precedence.

**Rationale**
- **Ridge** processes create new plate boundary; they should suppress subduction/transform in the same lattice cell.
- **Subduction** dominates over **transform** where convergence normal component is non‑trivial.

**Algorithm**
1. Build candidacy rasters **R**, **S**, **T** with confidence scores in [0,1].
2. Apply precedence filter to produce **M_final**:
   
   **M = pick(R>τ_R ? R : 0)**
   
   **M = pick(S>τ_S & M==0 ? S : M)**
   
   **M = pick(T>τ_T & M==0 ? T : M)**
3. Optionally dilate **M** by ≤1 cell for continuity (careful with narrow straits).

**Code sketch**
```ts
function resolveExclusivity(R,S,T, tauR,tauS,tauT){
  const M = zerosLike(R);
  forEachCell(i => {
    if (R[i] > tauR) { M[i] = RIDGE; return; }
    if (S[i] > tauS) { M[i] = SUBD;  return; }
    if (T[i] > tauT) { M[i] = XFORM; return; }
  });
  return M;
}
```

---

## SB‑01 — Subduction obliquity weighting

**Why**: Discount highly oblique segments to avoid spurious arcs/back‑arc edits and over‑accretion.

**Theory**
- Let **θ = arccos( (v·n̂)/\|v\| )**. Weight convergence by **w_obl = sin(θ)^p** (use **p=1.5** as proposed).
- Effective normal convergence used for arc/forearc/back‑arc rates:
  
  **v_eff = max(0, (v·n̂)) · w_obl**.

**Boundaries**
- θ≈0° (orthogonal) ⇒ w≈0 (no obliquity penalty). θ≈90° ⇒ w≈1.

---

## SB‑02 — Convergence gate (τ_conv) — raised threshold

**Why**: You raised the gate to **τ_conv ≥ 0.020 m/yr** to trim weak/oblique segments. Keep this as a *hard gate* before any subduction edits.

**Check**
- Use **v_eff** (from SB‑01). If **v_eff < τ_conv** then **skip subduction edits** for that cell this step.

---

## SB‑03 — Soft caps on arc/back‑arc per‑Myr edits

**Why**: Prevent runaway edits; keep rates within plausible geologic envelopes.

**Policy**
- **Arc**: ≤ **8 m/Myr**. **Back‑arc**: ≤ **6 m/Myr** (values from your guardrails).
- Apply the cap *after* obliquity weighting and τ_conv gating.

**Code**
```ts
arcRate = clamp(arcRateRaw, -ARC_CAP, ARC_CAP); // m/Myr
backArcRate = clamp(backArcRateRaw, -BACKARC_CAP, BACKARC_CAP);
```

---

## SB‑04 — Accretion & forearc mass balance sanity

**Intent**: Keep sedimentary/tectonic accretion increments small, tied to **v_eff** and taper geometry.

**Heuristic**
- d(th_c)_acc ∝ v_eff · f(taper, sediment)
- Keep cumulative **dC (mass)** within 0–0.003 per Myr (your logs show ≤0.0025), unless docked terranes are explicitly modeled.

---

## RF‑01 — Oceanization kinetics (continental → oceanic)

**Why**: You increased **k_c_oceanize** from 0.01 to **0.03 /Myr** to move more cores fully oceanic under sustained rifting.

**Model**
- Track an oceanization state **χ ∈ [0,1]**. Update as:
  
  **dχ/dt = k_c_oceanize · A_rift · (1 − χ)**,
  
  where **A_rift** is a rift activity factor in [0,1] driven by extensional strain rate.
- Switch material when **χ ≥ χ* (≈0.9)**.

**Note**
- Cap **k_c_oceanize·dt ≤ 0.2** for stability (dimensionless “reaction CFL”).

---

## RF‑02 — Rifting thinning from extensional strain rate

**Theory**
- If **ε̇** is the 2‑D extensional strain rate along the rift axis, approximate lithospheric thinning:
  
  **d(th_c)/dt ≈ −β · ε̇ · th_c**, with **β ∈ [0.7,1.0]**.
- Convert from plate kinematics: **ε̇ ≈ (Δv_parallel / L)** within the rift band.

**Guard**
- Keep **|d(th_c)/dt| ≤ 150 m/Myr** (current cap). Consider making the cap **derived from ε̇** instead of absolute.

---

## RF‑03 — Post‑rift thermal subsidence (oceanic age–depth)

**Background**
- Half‑space cooling gives **depth(age) ≈ d₀ + k√age** up to ~70 Ma.
- Plate model asymptote beyond ~70–80 Ma.

**Use**
- For diagnostics, compute **residual**: **r = z_model − z_ref(age)** and track mean/std/max.
- Targets: mean |r| ≲ 300–500 m for plausible maps; handle outliers where age resets.

---

## TR‑01 — Transform classification & seeding (releasing vs restraining)

**Theory**
- Decompose relative velocity at boundary into tangential **t̂** and normal **n̂**.
- **Releasing (pull‑apart)** if **(v·n̂) > −τ_tiny** and \|t\|>τ_t; **Restraining** if **(v·n̂) < −τ_tiny** with significant \|t\|.

**Implementation**
- Use your diagnostic: report **|t| min/mean/max** and **|n|**. Seed transforms where \|t\| exceeds percentile‑based threshold and **|n|** is small.

---

## TR‑02 — Suppress normal leakage in transform bands

**Why**: Prevent unintended opening/closing while applying shear.

**Method**
- In transform band pass, set **v_n = 0** (or blend to zero by Gaussian weight to avoid discontinuity) before shear advection.

---

## BZ‑01 — Buoyancy trust region (per‑step) & isostasy basis

**Current**: You log **dZ_step max = 200 m** and **amp ≈ 1.95 km**.

**Theory (Airy isostasy)**
- Vertical adjustment from crustal thickness change **Δh_c**:
  
  **Δz ≈ (ρ_m − ρ_c)/ρ_m · Δh_c**.
- For mantle lithosphere changes, include a term with **ρ_l** accordingly.

**Trust region**
- Apply **|Δz_step| ≤ Z_cap_step** (e.g., 200 m) and integrate with sub‑steps if needed so that **Σ|Δz_step|** respects stability and avoids ringing.

**Boundary‑prox scaling** (see BZ‑02) reduces buoyancy near active trenches/transforms.

---

## BZ‑02 — Boundary‑prox buoyancy down‑weighting

**Why**: Avoid double‑counting deformation that is already prescribed by boundary processes.

**Weighting**
- Let **d(x)** be distance to nearest active trench/transform. Define a taper:
  
  **w_b(d) = clamp((d − d₀)/(d₁ − d₀), 0, 1)**
  
  with **d₀ ≈ 0 km**, **d₁ ≈ 150–300 km**.
- Multiply buoyancy uplift/subsidence by **w_b(d)**.

---

## DX‑01 — Step diagnostics

**Add**
- Percent of boundary length by class (ridge/subduction/transform).
- Counts of newly oceanized cells, and totals.
- Top‑k CFL offenders per class and recommended **N_sub**.
- Hypsometry + land fraction.

**Example**
```txt
[diag] boundaries: ridge=27%, subd=41%, transform=32% | oceanized new=14 total=932
```

---

## RN‑01 — Overlay projection into raster rect

**Why**: You’ve fixed the offset by projecting overlays into the actual raster rect. Keep this invariant: **UI overlays sample the *same* raster transform as the simulation**.

**Checklist**
- One source of truth for the world→raster transform.
- Avoid per‑layer bespoke resampling.

---

## CH‑01 — Dimensionless stepping guards

**Principle**: express rate caps as *dimensionless* constraints times a scale.

- Shear CFL: **vt·dt / w_half ≤ 0.5**.
- Reaction CFL (oceanization): **k·dt ≤ 0.2**.
- Geomorphic caps per‑Myr: keep under **10 m/Myr** unless justified.

---

## CH‑02 — Deterministic pass order & reproducibility

**Guidelines**
- Pure functions for each pass. No hidden frame state.
- Order: classify→resolve exclusivity→compute sub‑steps per class→apply passes (ridge/subduction/transform)→secondary (accretion, arcs, back‑arc)→buoyancy→erosion/flexure.
- Seed RNG with the run seed; never call RNG inside geometry kernels unless seeded sub‑step locally.

---

## CH‑03 — Constants, tuning knobs & logging

**Do**
- Collect tunables in a single typed object, e.g., `ModelParams`.
- Unit‑suffixed names: `tauConv_m_per_yr`, `wHalfTransform_m`.
- Log **pre** and **post** values when clamping (you already print caps—great).

---

## VL‑01 — Sanity targets & calibration numbers

- **Plate speeds**: 20–100 mm/yr ⇒ 0.02–0.10 m/yr typical. If you see vt≈0.7 m/yr, check scaling.
- **Hypsometry**: land fraction ~25–35% (Earth today ~29%); for stylized worlds we accept 5–40% by design.
- **Rift thinning**: ≤100–200 m/Myr are within plausible first‑order bulk lithosphere changes.
- **Orogen thickening**: 200–400 m/Myr is conservative at crustal thickness scale (vertical surface uplift often 1–10 mm/yr but not 1:1 with thickening).
- **Age–depth residual**: mean |r| < ~400 m, σ < ~1 km, check spikes > 6–8 km.

---

## QA‑01 — Minimal regression tests (CPU‑only)

- **Determinism**: same seed → identical raster hashes.
- **CFL**: assert no step reports CFL > 0.5 after sub‑stepping.
- **Exclusivity**: no cell has two primary classes in the same step.
- **Oceanization monotonicity**: χ is non‑decreasing under rifting; doesn’t change otherwise.

**Sketch**
```ts
it('cfl bounded', () => {
  const log = runWorld({steps:200});
  expect(max(log.cfl.values)).toBeLessThanOrEqual(0.5 + 1e-6);
});
```

---

## Appendix A — Useful formulas

- **Obliquity**: `theta = Math.acos( dot(v, nHat) / (norm(v)+eps) )` then `w = Math.pow(Math.sin(theta), 1.5)`.
- **Distance taper**: `w = clamp01((d - d0) / (d1 - d0))`.
- **Half‑space cooling (diagnostic)**: `depth = d0 + k * Math.sqrt(age_Ma)`; choose `d0≈2600 m`, `k≈350 m/√Ma`.
- **Isostasy (crust)**: `dz = ((rho_m - rho_c)/rho_m) * d_th_c`.
- **Sub‑step count**: `N = ceil( (v*dt) / (CFLmax * w_half) )`.

---

## Appendix B — Implementation notes for GPU/CPU raster coupling

- Keep one **face→bary→lattice** mapping function; test with a checkerboard.
- Export ROI probes deterministically (you already have `rollover probe CSV`).
- Avoid read‑modify‑write hazards across passes: use double‑buffered rasters or accumulate deltas then commit.

---

## Appendix C — Logging lines (reference)

- ` [cfl] class: vt*dt/w_half = X > 0.5 (vt=..., dt=..., w_half=...)`
- ` [substep] class: increasing substeps A → B (CFL Y)`
- ` [cap] {process}: rate capped (raw=..., cap=...)`
- ` [diag] age-depth residual: mean=..., std=..., max|r|=...`
- ` [diag] hyps: land=...%`
- ` [transforms.diag] from_boundaries=..., passing=..., |t| min/mean/max=..., |n| min/mean/max=...`

---

**End of document.**

