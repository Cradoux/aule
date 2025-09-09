Awesome—thanks for the zip. I read through the engine crate (boundary finding, subduction, rifting, transforms, accretion/orogeny, flexure/isostasy, pipeline + diagnostics). Below is a tight but thorough review from two angles: (A) simulation fundamentals, so we can avoid hacks; (B) code quality/architecture so the model is easier to reason about and tune.

# Simulation review (what the logs + code are telling us)

**1) Rates are routinely orders-of-magnitude above physical and then hard-clamped.**
You’re seeing raw thickening 300–600 m/Myr (even 600+) and rift thinning \~700–950 m/Myr, then clamping at 150–300 m/Myr every step. That means the *rate laws* are the real problem; the caps are masking it. Same for transforms: CFL blows up (vt\*dt/w ≳ 2–3) and you react with 32× substeps. This is numerically stable, but it’s compensating for rate formulas that are too aggressive for the chosen time step and widths.

**Action:** derive each tectonic rate from kinematics + characteristic width rather than “free” coefficients.

* **Transforms:** set shear strain rate $\dot\gamma \approx v_t/W_t$. Pick $W_t$ that scales with lithospheric thickness/age or a calibrated constant (e.g., 20–50 km). Then compute displacement per step as $v_t\,\Delta t$ and ensure CFL \~ 0.2–0.4 by design. Your current hardcoded $W_{half}=10\,\text{km}$ is likely too small for your dt and vt.
* **Convergence/orogeny:** let thickening rate follow $\dot{h}_c=\alpha\,v_n\,\sin^\beta\theta/W_c$, with $W_c$ a compressive belt width (100–300 km typical), $\alpha\in[0.3,1]$, $\beta\simeq 1\!-\!2$. This yields realistic 1–10 mm/yr vertical rates without post-hoc clamps.
* **Rifting:** tie thinning to horizontal extension: $\dot{h}_c=-(v_{rift}/W_r) \cdot h_c$, or implement a McKenzie‐style β evolution with a relaxation timescale (e.g., $d\beta/dt=(\beta_\text{target}-\beta)/\tau$, $\tau=5\!-\!15$ Myr). That will naturally cap rates and avoid huge negative spikes you’re clamping now.
* **Oceanization:** your new $k_{c,\text{oceanize}}=0.03/\text{Myr}$ aggressively drives many cores fully oceanic; with the above rift law you can lower this back toward 0.01–0.02 and let kinematics do the work.

**2) Obliquity and exclusivity should be first-class, not afterthoughts.**
You’re already gating subduction with τconv and you mentioned a sin(θ)^1.5 weight. Do it—and also enforce **mutually exclusive process masks per cell/edge** each step: ridge → subduction → transform. That single change will eliminate a lot of “stacked” edits and reduce your need for caps.

**3) Buoyancy/isostasy trust region is a symptom.**
The step clamp $dZ_{\text{step}}\le 200$ m and “amp \~1957 m” saturating suggest you’re repeatedly overshooting the vertical equilibration. That usually means (a) inconsistent densities or (b) too-large dt relative to your relaxation model.

* Make the **isostatic update implicit** (solve $Z^{k+1}$ from the load relation directly) or add a proper **line search**/backtracking so you don’t need a blunt 200 m cap.
* Down-weight buoyancy *inside active trench/transform bands* (what you proposed) but do it as a physical mask: e.g., set $W_b$ (isostatic support length) shorter under strong tectonic forcing, rather than a pure multiplier on ΔZ.

**4) Land fraction staying \~6%**
That’s the rift/oceanization dominance + frequent thickening caps. Calibrate the above rates and you should organically climb toward \~25–35% without artificial sea-level tweaks.

**5) Diagnostics to make the rates honest**
Great that you log τconv, |t|, |n|, etc. Add:

* boundary length % per class (you listed),
* fraction of cells hitting each cap,
* global mass/volume budgets (net crustal volume change vs addition/removal by arcs/back-arc/subduction),
* mean/max $|d h_c|$, $|dZ|$ per process *before* and *after* any limiter.

These will quickly tell you whether caps are doing heavy lifting (they are).

# Concrete modeling tweaks (low-risk, high-impact)

1. **Exclusivity masks (do it now).**
   Bitflags per cell/edge with a fixed priority: Ridge → Subduction → Transform. Apply once at the start of the step and consume flags as processes run.

2. **Subduction obliquity weighting.**
   Multiply your effective convergence by $\sin^\beta\theta$ with $\beta ≈ 1.5$. You already compute $\theta$; wire this into both interface fluxes and any arc/back-arc budgets. This will trim those oblique segments that still sneak through your τconv gate.

3. **Boundary-prox buoyancy scaling.**
   Define a distance‐to‐boundary field $d(x)$ and reduce isostatic response inside $d<W_b$ with a smooth kernel, e.g. $w_b(d)=\tfrac{1}{2}[1+\tanh((d-d_0)/s)]$. Use $d_0≈50$ km, $s≈20$ km. This is a *physical* de-emphasis where tectonic forcing dominates short-term vertical motion.

4. **Unified CFL limiter** (replace per-module ad-hoc substeps):
   For any process that moves material a distance $u\,\Delta t$ across a width $W$, enforce $C=u\Delta t/W\le C_{max}$ by *scaling the update*, not only by increasing substeps. This keeps updates single-pass and predictable.

```rust
fn cfl_limited_update(raw: f64, width: f64, dt: f64, cfl_max: f64) -> f64 {
    // raw is the intended displacement (m) over dt
    let c = (raw.abs() / width).max(1e-12);
    let scale = (cfl_max / c).min(1.0);
    raw * scale
}
```

5. **Replace hard caps with smooth saturators** so derivatives don’t jump and you stop fighting the integrator:

```rust
#[inline]
fn soft_cap(x: f64, cap: f64) -> f64 {
    let a = x / cap;
    cap * a.tanh()        // odd, smooth, ≤ cap in magnitude
}
```

Use this once per process on the *rate*, not the accumulated Δ over dt.

# Code/architecture review (to make the physics easy to tune)

**Units, everywhere.**
You mix m, km, Myr in locals—correctly—but it’s easy to slip. Adopt `uom` (units-of-measure) or a thin newtype layer (Meters, Myr, MPerYr). That alone prevents half the “why is CFL exploding” bugs.

**Centralize knobs in a typed `SimParams`.**
Right now many magic numbers (0.025, 10\_000.0, 150.0 m/Myr, etc.) live in modules. Hoist them into a single struct with doc comments + units, and read-only copies per process for hot loops. Add a `tuning` profile (conservative/standard/aggressive) so experiments are reproducible.

**One scheduler, many processes.**
Pipeline order is clear, but substepping is scattered. Put a tiny “process runner” that:

* asks each process for its local CFL and suggested width,
* computes a global per-step scale (or substeps),
* applies the scaled update in a consistent order.
  This removes the repeated “increasing substeps 4 → 32 (CFL …)” spam and makes performance predictable.

**Exclusivity as bitflags.**
Represent boundary/process eligibility as `bitflags!` on per-edge/per-cell arrays. Resolve once per step. This removes ad-hoc `if ridge && subduction` checks sprinkled around.

**Error handling & logging.**
You’re mostly clean of `unwrap!`, which is good. Consider:

* `tracing` with spans per process/step; throttle high-rate logs (your transform diag repeats dozens of times due to substeps).
* Summaries once per Myr with means/stddevs + % capped.
* Guard rails: assert invariants (no negative thickness; |slope| < …; mass budgets within tolerance).

**Numerics.**
Use f64 for accumulations even if storage is f32. You do this in many places; keep it consistent. For spherical geometry, your `local_basis` has a good polar fallback—nice—just keep it the single source of truth (I saw a few local re-builds of tangents; funnel them through geo.rs).

**Testing.**
Add a few extremely cheap tests:

* Property tests: isostasy returns zero ΔZ when loads match reference; exclusivity never double-applies; CFL limiter never exceeds target.
* Snapshot tests for one step on a tiny grid (8×4): store boundary counts, land %, cap hit %.

# Priorities I’d hit next (in order)

1. **Implement exclusivity masks + obliquity weight** (fast, surgical, removes stacking).
2. **Replace per-module caps with the `soft_cap` and a *unified CFL limiter* based on physical widths** (stops the log spam and substep explosions).
3. **Re-derive rift/thicken laws from $v/W$** and calibrate $W_t, W_r, W_c$ so raw rates land in 10–200 mm/yr, not 900.
4. **Implicit or line-searched isostasy update** so the “200 m trust region” disappears.
5. **Boundary-prox buoyancy weighting** (your idea—keep it smooth and width-based).
6. **Diagnostics you listed** (+ mass/volume budgets).

If you want, I can sketch exact formulas for the three widths and give starter parameter values so that, with your current dt, CFL stays ≤0.3 without substeps and land fraction trends toward \~30% after \~60 Myr—*without* any hard caps.
