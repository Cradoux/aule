## Continental vs Oceanic: Fields, Updates, and Classification

This document summarizes how the engine distinguishes and evolves continental vs oceanic domains.

Key fields (World): `c` (continental fraction 0..1), `th_c_m` (continental crust thickness, m), `age_myr`, `area_m2`, and `plates.plate_id`.

- Initialization: `continent::build_continents` builds a deterministic template. When continents are enabled in the pipeline, `c` is seeded from that template (0..1) and `th_c_m` initialized to a uniform base (e.g., 2500 m).

- Advection: `continent::advect_c_thc` semi-Lagrangian advection on the grid moves `c` and `th_c_m` with rigid plate motion. The advection now clamps `c` to [0,1] and caps per-step `Δth_c` to ±200 m to avoid unit mistakes.

- Modification Operators:
  - Rifting (`rifting.rs`): thins crust and reduces `c` in extension bands (oceanization). Uses plate kinematics and boundary masks.
  - Accretion (`accretion.rs`) and Orogeny (`orogeny.rs`): increase `th_c_m` and (often) `c` near arcs or collision zones. Per-operator caps and unit sanity logs are in place.
  - Surface processes (`surface.rs`): operate on topography and sediment. `c` is currently not altered here (future: weathering-derived changes could be added).

- Buoyancy Coupling: `continent::apply_uplift_from_c_thc` adjusts `depth_m` via isostatic-like coupling using `c` and `th_c_m` anomalies relative to a continental reference (35 km), scaled by buoyancy factor `(ρ_m - ρ_c)/ρ_m`.

- Gates/Thresholds: Several operators use `c_min` thresholds for continental behavior, e.g., subduction/orogeny commonly use `C ≥ 0.60` to gate continental responses. Transforms generally ignore `c` (shear-dominated).

- Missing or Addressed Advection: Earlier nearest-neighbor swaps for `th_c_m` produced km-scale deltas; replaced by 1-ring weighted interpolation with a ±200 m per-step cap. `c` is explicitly clamped to [0,1] after advection.

### Classification Rule for Plate Type Overlay

- Per-plate classification uses area-weighted mean C̄ over the plate:
  - Continental if `C̄ ≥ 0.60`
  - Oceanic if `C̄ ≤ 0.20`
  - Mixed otherwise

- A per-cell mode colors by `c` (0..1) with a continuous colormap.

### Code Pointers

- Fields: `engine/src/world.rs` (`c`, `th_c_m`, `area_m2`, `plates.plate_id`).
- Advection: `engine/src/continent.rs::advect_c_thc`.
- Rifting: `engine/src/rifting.rs` (thinning/oceanization).
- Accretion/Orogeny: `engine/src/accretion.rs`, `engine/src/orogeny.rs`.
- Buoyancy: `engine/src/continent.rs::apply_uplift_from_c_thc`.
- Pipeline: `engine/src/pipeline.rs` (rigid motion, advection, operators, clamping).
- Viewer overlays: `viewer/src/overlay.rs` (Plate ID points, Plate Type overlay), UI in `viewer/src/main.rs`.


