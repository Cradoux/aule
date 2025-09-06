//! World state container and constructors (T-400).
//!
//! Conventions and scope:
//! - `depth_m` uses positive-down bathymetry; land elevation is `-depth_m`.
//! - Sea-level regulation (isostasy/eustasy) is modeled as a uniform offset added to `depth_m`.
//! - Age→depth sets an oceanic baseline; tectonic edits (subduction, transforms, rifting, orogeny)
//!   and continent uplift then modify `depth_m` on top.
//! - Flexure is approximated by a Winkler-like response in CPU MVP; GPU solver lives elsewhere.
//!
//! Important integration notes:
//! - There are two stepping paths: this file's `step_once` (CPU pipeline that applies isostasy by
//!   writing back into `depth_m`) and `engine::pipeline::step_full` (viewer pipeline) which solves
//!   an `eta` but does not modify `depth_m`. Callers that use the pipeline must render with
//!   `elev = -depth - eta` to get consistent coastlines.

use crate::{
    boundaries::Boundaries, continent, flexure_loads, grid::Grid, isostasy, plates::Plates, ridge,
    subduction, transforms,
};
use std::time::Instant;

/// Simulation clock information.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Clock {
    /// Current simulation time in Myr.
    pub t_myr: f64,
    /// Step index (starts at 0, increments per step).
    pub step_idx: u64,
}

/// The complete world state required for stepping.
pub struct World {
    /// Geodesic grid definition.
    pub grid: Grid,
    /// Plate assignment and Euler poles.
    pub plates: Plates,
    /// Classified boundaries for current state.
    pub boundaries: Boundaries,
    /// Age field in Myr, length = grid.cells.
    pub age_myr: Vec<f32>,
    /// Staging buffer for oceanic age during rigid advection
    pub age_stage: Vec<f32>,
    /// Bathymetry depth in meters (+ down), length = grid.cells.
    pub depth_m: Vec<f32>,
    /// Staged depth for end-of-step atomic publish (same length as depth_m)
    pub depth_stage_m: Vec<f32>,
    /// Continental fraction (0..1). Seeded when continents are enabled.
    pub c: Vec<f32>,
    /// Staging buffer for continental fraction during rigid advection
    pub c_stage: Vec<f32>,
    /// Continental crust thickness in meters (≥0). Seeded alongside `C`.
    pub th_c_m: Vec<f32>,
    /// Staging buffer for continental crust thickness during rigid advection
    pub th_c_stage: Vec<f32>,
    /// Per-cell velocities (east,north) in m/yr.
    pub v_en: Vec<[f32; 2]>,
    /// Simulation clock.
    pub clock: Clock,
    /// Reference sea-level volume and area captured after baseline age→depth
    pub sea_level_ref: Option<SeaLevelRef>,
    /// Precomputed per-cell areas in m^2 on the sphere (4πR^2 scaled)
    pub area_m2: Vec<f32>,
    /// Last flexure residual ratio (r_out / max(1, r_in)) if flexure applied this step.
    pub last_flex_residual: f32,
    /// Cumulative sediment thickness in meters (≥0), updated by surface processes.
    pub sediment_m: Vec<f32>,
    /// Elastic thickness field Te (meters), used to parameterize flexure response.
    pub te_m: Vec<f32>,
    /// Persistent buoyancy delta amplitude (meters), composed into depth each step
    pub delta_buoy_m: Vec<f32>,
    /// Last surface-processes stats (if applied this step).
    pub last_surface_stats: Option<crate::surface::SurfaceStats>,
    /// Continents change epoch counter (bump when C/th_c change applied).
    pub epoch_continents: u64,
    /// Last epoch when sea-level was re-baselined (to debounce auto mode).
    pub last_rebaseline_epoch: u64,
    /// Sea state (offset from geoid). Positive raises sea level (more ocean).
    pub sea: SeaState,
    /// Last-step crustal mass (kg) for budget deltas.
    pub last_mass_cont_kg: f64,
    /// Last-step sediment mass (kg) for budget deltas.
    pub last_mass_sed_kg: f64,
    /// Last-step ocean volume (m^3) for budget deltas.
    pub last_ocean_vol_m3: f64,
    /// Persistent scratch buffers for pipeline to reduce transient allocations.
    pub scratch: WorldScratch,
}
/// Reusable scratch space to avoid per-step allocations (sizes match grid.cells).
pub struct WorldScratch {
    /// Scratch float buffer A, length equals `grid.cells`. Caller defines semantics per stage.
    pub f32_a: Vec<f32>,
    /// Scratch float buffer B, length equals `grid.cells`. Caller defines semantics per stage.
    pub f32_b: Vec<f32>,
    /// Scratch float buffer C, length equals `grid.cells`. Caller defines semantics per stage.
    pub f32_c: Vec<f32>,
    /// Scratch float buffer D, length equals `grid.cells`. Caller defines semantics per stage.
    pub f32_d: Vec<f32>,
}

/// Simple-mode preset parameters for deterministic world generation.
#[derive(Clone, Copy, Debug)]
pub struct SimplePreset {
    /// Number of plates to seed deterministically
    pub plates: u32,
    /// Number of synthetic continental caps
    pub continents_n: u32,
    /// Target land fraction in [0,1]
    pub target_land_frac: f32,
}

/// Summary of the simple generation pass.
#[derive(Clone, Copy, Debug)]
pub struct SimpleReport {
    /// Area-weighted ocean fraction (0..1)
    pub ocean_frac: f32,
    /// Min/max of depth (m); depth>0 water, depth<0 land elevation = -depth
    pub depth_min_max: (f32, f32),
    /// Area-weighted land fraction (0..1)
    pub land_frac: f32,
}

/// Reference sea-level bookkeeping
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeaLevelRef {
    /// Target ocean volume at reference state (m^3)
    pub volume_m3: f64,
    /// Ocean area used at reference state (m^2)
    pub ocean_area_m2: f64,
}

/// Sea state parameters (viewer-controlled for Simple mode)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeaState {
    /// Sea-level offset η in meters relative to geoid. Positive is deeper water.
    pub eta_m: f32,
}

impl World {
    // removed old helpers in favor of sea_level module functions

    #[allow(dead_code)]
    /// Evaluate land fraction after rebaseline to a fixed target water volume.
    fn land_frac_after_rebaseline_to_volume(
        &self,
        base_depth: &[f32],
        cont_mask: &[f32],
        amp_m: f32,
        target_volume_m3: f64,
    ) -> (f32, f64) {
        let mut depth = base_depth.to_vec();
        for i in 0..depth.len() {
            depth[i] -= amp_m * cont_mask[i];
        }
        let off =
            isostasy::solve_offset_for_volume(&depth, &self.area_m2, target_volume_m3, 1e6, 64);
        let mut land_area = 0.0f64;
        let mut total_area = 0.0f64;
        for (d, a) in depth.iter().zip(self.area_m2.iter()) {
            let val = (*d as f64) + off;
            total_area += *a as f64;
            if val < 0.0 {
                land_area += *a as f64;
            }
        }
        let frac = if total_area > 0.0 { (land_area / total_area) as f32 } else { 0.0 };
        (frac, off)
    }

    #[allow(dead_code)]
    /// Bisection on amplitude while locking water volume to a fixed target.
    fn solve_amplitude_for_land_fraction_locked_volume(
        &mut self,
        cont_mask: &[f32],
        target_land: f32,
        target_volume_m3: f64,
        amp_hi_m: f32,
        tol_land: f32,
    ) -> f32 {
        let base_depth = self.depth_m.clone();
        let mut lo = 0.0f32;
        let mut hi = amp_hi_m;
        let mut _f_lo = {
            let (lf, _) = self.land_frac_after_rebaseline_to_volume(
                &base_depth,
                cont_mask,
                lo,
                target_volume_m3,
            );
            lf - target_land
        };
        let mut _f_hi = {
            let (lf, _) = self.land_frac_after_rebaseline_to_volume(
                &base_depth,
                cont_mask,
                hi,
                target_volume_m3,
            );
            lf - target_land
        };
        if _f_lo.abs() <= tol_land {
            // Already matching at zero amplitude
            self.depth_m = base_depth;
            return lo;
        }
        if _f_hi <= 0.0 {
            // Even huge uplift cannot reach target: clamp to hi
            let (_lf, off) = self.land_frac_after_rebaseline_to_volume(
                &base_depth,
                cont_mask,
                hi,
                target_volume_m3,
            );
            self.depth_m.clone_from(&base_depth);
            for (i, &cm) in cont_mask.iter().enumerate().take(self.depth_m.len()) {
                self.depth_m[i] -= hi * cm;
            }
            self.sea.eta_m = off as f32;
            return hi;
        }
        for _ in 0..48 {
            let mid = 0.5 * (lo + hi);
            let (lf, off) = self.land_frac_after_rebaseline_to_volume(
                &base_depth,
                cont_mask,
                mid,
                target_volume_m3,
            );
            let f_mid = lf - target_land;
            if f_mid.abs() <= tol_land {
                self.depth_m.clone_from(&base_depth);
                for (i, &cm) in cont_mask.iter().enumerate().take(self.depth_m.len()) {
                    self.depth_m[i] -= mid * cm;
                }
                self.sea.eta_m = off as f32;
                return mid;
            }
            if f_mid > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let mid = 0.5 * (lo + hi);
        let (_lf, off) = self.land_frac_after_rebaseline_to_volume(
            &base_depth,
            cont_mask,
            mid,
            target_volume_m3,
        );
        self.depth_m.clone_from(&base_depth);
        for (i, &cm) in cont_mask.iter().enumerate().take(self.depth_m.len()) {
            self.depth_m[i] -= mid * cm;
        }
        self.sea.eta_m = off as f32;
        mid
    }
    /// Construct a world with deterministic plates and zeroed ages/depths.
    pub fn new(f: u32, num_plates: u32, seed: u64) -> Self {
        let grid = Grid::new(f);
        // Precompute local bases/edge lengths once at world creation
        let mut grid = grid;
        grid.precompute_local_bases_and_lengths();
        let plates = Plates::new(&grid, num_plates, seed);
        let v_en = plates.vel_en.clone();
        let boundaries =
            crate::boundaries::Boundaries::classify(&grid, &plates.plate_id, &v_en, 0.005);
        let age_myr = vec![0.0f32; grid.cells];
        let age_stage = vec![0.0f32; grid.cells];
        let depth_m = vec![0.0f32; grid.cells];
        let depth_stage_m = vec![0.0f32; grid.cells];
        let c = vec![0.0f32; grid.cells];
        let c_stage = vec![0.0f32; grid.cells];
        let th_c_m = vec![35_000.0f32; grid.cells];
        let th_c_stage = vec![35_000.0f32; grid.cells];
        let sediment_m = vec![0.0f32; grid.cells];
        let cells = grid.cells;
        let clock = Clock { t_myr: 0.0, step_idx: 0 };
        // Precompute area in m^2 using Earth radius
        const R_EARTH_M: f64 = 6_371_000.0;
        let scale = 4.0 * std::f64::consts::PI * R_EARTH_M * R_EARTH_M;
        let mut area_m2: Vec<f32> = Vec::with_capacity(grid.cells);
        for &a in &grid.area {
            area_m2.push((a as f64 * scale) as f32);
        }
        Self {
            grid,
            plates,
            boundaries,
            age_myr,
            age_stage,
            depth_m,
            depth_stage_m,
            c,
            c_stage,
            th_c_m,
            th_c_stage,
            v_en,
            clock,
            sea_level_ref: None,
            area_m2,
            last_flex_residual: 0.0,
            sediment_m,
            te_m: vec![25_000.0; cells],
            delta_buoy_m: vec![0.0; cells],
            last_surface_stats: None,
            epoch_continents: 0,
            last_rebaseline_epoch: 0,
            sea: SeaState { eta_m: 0.0 },
            last_mass_cont_kg: 0.0,
            last_mass_sed_kg: 0.0,
            last_ocean_vol_m3: 0.0,
            scratch: WorldScratch {
                f32_a: vec![0.0; cells],
                f32_b: vec![0.0; cells],
                f32_c: vec![0.0; cells],
                f32_d: vec![0.0; cells],
            },
        }
    }

    /// Deterministic end-to-end generation for Simple mode.
    pub fn generate_simple(&mut self, preset: &SimplePreset, seed: u64) -> SimpleReport {
        // 1) Reset/normalize world state for current F & plates
        self.plates = crate::plates::Plates::new(&self.grid, preset.plates, seed);
        self.v_en.clone_from(&self.plates.vel_en);
        self.age_myr.fill(0.0);
        self.age_stage.clone_from(&self.age_myr);
        self.depth_m.fill(0.0);
        self.c.fill(0.0);
        self.c_stage.clone_from(&self.c);
        self.th_c_m.fill(35_000.0);
        self.th_c_stage.clone_from(&self.th_c_m);
        // Enforce continental thickness bounds after init
        for t in &mut self.th_c_m {
            *t = t.clamp(25_000.0, 65_000.0);
        }
        self.sediment_m.fill(0.0);
        self.te_m.fill(25_000.0);
        self.delta_buoy_m.resize(self.grid.cells, 0.0);
        self.delta_buoy_m.fill(0.0);
        self.clock = Clock { t_myr: 0.0, step_idx: 0 };
        self.sea_level_ref = None;
        self.last_mass_cont_kg = 0.0;
        self.last_mass_sed_kg = 0.0;
        self.last_ocean_vol_m3 = 0.0;
        self.boundaries =
            Boundaries::classify(&self.grid, &self.plates.plate_id, &self.v_en, 0.005);

        // Derive per-plate kind from initial continental mask (C field)
        self.plates.kind = crate::plates::derive_kinds(
            &self.grid,
            &self.plates.plate_id,
            &self.c,
            0.35, // plate is Continental if >35% of its area has C>0.5
            0.5,
        );

        // 2) Initialize ridge births then baseline bathymetry from age (ocean depths)
        {
            let mut ages_tmp = self.age_myr.clone();
            let _ridge_stats = crate::ridge::apply_ridge(
                &self.grid,
                &self.boundaries,
                &mut ages_tmp,
                crate::ridge::RidgeParams { fringe_age_myr: 0.0 },
            );
            self.age_myr = ages_tmp;
        }
        let n_cells = self.grid.cells;
        for i in 0..n_cells {
            let mut d = crate::age::depth_from_age_plate(
                self.age_myr[i] as f64,
                2600.0,
                self.clock.t_myr,
                6000.0,
                1.0e-6,
            ) as f32;
            if !d.is_finite() {
                d = 6000.0;
            }
            self.depth_m[i] = d.clamp(0.0, 6000.0);
        }
        // Keep stage in sync after initialization
        self.depth_stage_m.clone_from(&self.depth_m);
        // Capture reference ocean volume for constant-volume isostasy using elevation definition
        let ref_sea = crate::isostasy::compute_ref(&self.depth_m, &self.area_m2, self.sea.eta_m);
        self.sea_level_ref = Some(ref_sea);

        // 3) Continents: build template.
        // If continents_n == 0, treat as a "supercontinent" preset and synthesize a single
        // contiguous ribbon with several lobes for early Wilson-cycle scenarios.
        let mut tpl: Vec<f32> = if preset.continents_n == 0 {
            crate::continent::build_supercontinent_template(
                &self.grid, seed, 5,      // lobes
                2200.0, // lobe radius (km)
                600.0,  // falloff (km)
            )
        } else {
            let cp = crate::continent::ContinentParams {
                seed,
                n_continents: preset.continents_n,
                mean_radius_km: 2200.0,
                falloff_km: 600.0,
                plateau_uplift_m: 1.0,
                target_land_fraction: None,
            };
            crate::continent::build_continents(&self.grid, cp).uplift_template_m
        };
        // Trim tiny template tails to avoid near-global uplift when solving amplitude
        let mut vmax = 0.0f32;
        for &v in &tpl {
            if v > vmax {
                vmax = v;
            }
        }
        if vmax > 0.0 {
            let thr = 0.05f32 * vmax;
            for v in &mut tpl {
                if *v < thr {
                    *v = 0.0;
                }
            }
        }

        // 4) Compute offset that would hit target ocean area on the baseline, and evaluate its ocean fraction
        let target_ocean_frac = (1.0 - preset.target_land_frac as f64).clamp(0.0, 1.0);
        let off_area = crate::sea_level::solve_offset_for_ocean_area_fraction(
            &self.depth_m,
            &self.area_m2,
            target_ocean_frac as f32,
        );
        let _ocean_at_off = {
            let mut ocean_area = 0.0f64;
            let mut tot = 0.0f64;
            for (d, a) in self.depth_m.iter().zip(self.area_m2.iter()) {
                if (*d as f64 + off_area) > 0.0 {
                    ocean_area += *a as f64;
                }
                tot += *a as f64;
            }
            if tot > 0.0 {
                ocean_area / tot
            } else {
                0.0
            }
        };

        // Optional inherited belts for supercontinent: shallow, wide bands
        if preset.continents_n == 0 {
            // Derive a deterministic seed for belts from stable world properties
            let seed_belts = (self.grid.frequency as u64)
                ^ ((self.plates.pole_axis.len() as u64) << 16)
                ^ (self.clock.step_idx);
            let belts = crate::continent::build_supercontinent_belts(
                seed_belts,
                crate::continent::BeltParams::default(),
            );
            let imprint = crate::continent::imprint_orogenic_belts(&self.grid, &belts);
            for (d, di) in self.depth_m.iter_mut().zip(imprint.into_iter()) {
                *d += di;
            }
        }

        // 5) Use principled amplitude solver only (no fixed/quantile/fallback paths)
        let (amp_m, off_m) = crate::isostasy::solve_amplitude_for_land_fraction(
            &tpl,
            &self.depth_m,
            &self.area_m2,
            preset.target_land_frac,
            0.0,
            6000.0,
            2e-2,
            1e-4,
            48,
        );
        for (i, &cm) in tpl.iter().enumerate().take(self.grid.cells) {
            self.depth_m[i] -= (amp_m as f32) * cm;
        }
        self.sea.eta_m = -(off_m as f32);

        // 6) Initialize continental fields so continents make up the target land mass
        // Build an area-targeted binary C from the template, matching preset.target_land_frac
        {
            let total_area: f64 = self.area_m2.iter().map(|&a| a as f64).sum();
            let target_area: f64 = (preset.target_land_frac as f64) * total_area;
            // Create sortable list of (tpl, area, idx)
            let mut cells: Vec<(f32, f32, usize)> = tpl
                .iter()
                .zip(self.area_m2.iter())
                .enumerate()
                .map(|(i, (t, a))| (*t, *a, i))
                .collect();
            // Sort by template strength descending
            cells.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            // Fill C to hit target area; allow fractional fill on the boundary cell
            let mut acc: f64 = 0.0;
            self.c.fill(0.0);
            for (val, a, idx) in cells {
                if val <= 0.0 {
                    break;
                }
                if acc >= target_area {
                    break;
                }
                let a64 = a as f64;
                let remaining = (target_area - acc).max(0.0);
                if remaining >= a64 {
                    self.c[idx] = 1.0;
                    acc += a64;
                } else {
                    self.c[idx] = (remaining / a64) as f32; // partial fill to hit target exactly
                    acc = target_area;
                }
            }
            // Seed a uniform continental thickness for C>0
            for (i, c) in self.c.iter().enumerate().take(self.grid.cells) {
                self.th_c_m[i] = if *c > 0.0 { 40_000.0 } else { 0.0 };
            }
            // Re-derive plate kinds now that C exists
            self.plates.kind =
                crate::plates::derive_kinds(&self.grid, &self.plates.plate_id, &self.c, 0.35, 0.5);
        }

        // 7) Skip subduction/flexure in Simple; keep land fraction stable

        // 8) No robustness fallback: keep the solver result deterministically

        // 9) Final: keep η managed externally (viewer) in Simple mode; do not modify depths here
        // Capture reference ocean state AFTER continents and final η using elevation definition
        let ref_after = crate::isostasy::compute_ref(&self.depth_m, &self.area_m2, self.sea.eta_m);
        self.sea_level_ref = Some(ref_after);

        // Report stats
        let mut dmin = f32::INFINITY;
        let mut dmax = f32::NEG_INFINITY;
        let mut dsum = 0.0f64;
        for &d in &self.depth_m {
            if d.is_finite() {
                if d < dmin {
                    dmin = d;
                }
                if d > dmax {
                    dmax = d;
                }
                dsum += d as f64;
            }
        }
        let mut ocean_area = 0.0f64;
        let mut total_area = 0.0f64;
        for (d, a) in self.depth_m.iter().zip(self.area_m2.iter()) {
            total_area += *a as f64;
            if (*d as f64 + self.sea.eta_m as f64) > 0.0 {
                ocean_area += *a as f64;
            }
        }
        let ocean_frac = if total_area > 0.0 { ocean_area / total_area } else { 0.0 };
        let land_frac = 1.0 - ocean_frac;
        println!(
            "[simple] land_target={:.0}% land={:.1}% depth[min/mean/max]={:.0}/{:.0}/{:.0} m",
            preset.target_land_frac as f64 * 100.0,
            land_frac * 100.0,
            dmin,
            (dsum / (self.depth_m.len().max(1) as f64)),
            dmax
        );
        SimpleReport {
            ocean_frac: ocean_frac as f32,
            depth_min_max: (dmin, dmax),
            land_frac: land_frac as f32,
        }
    }
}

/// Parameters that control a single evolution step.
#[derive(Clone, Copy, Debug)]
pub struct StepParams {
    /// Time step in Myr
    pub dt_myr: f64,
    /// Apply elastic flexure response to current loads
    pub do_flexure: bool,
    /// Adjust global sea level to maintain reference ocean volume
    pub do_isostasy: bool,
    /// Apply transform pull-apart/restraining bands
    pub do_transforms: bool,
    /// Apply subduction trench/arc/backarc edits
    pub do_subduction: bool,
    /// Advect C/th_c and apply continental uplift to depth
    pub do_continents: bool,
    /// Reset age along divergent boundaries (ridge births)
    pub do_ridge_birth: bool,
    /// If true, auto re-baseline sea level after continents change.
    pub auto_rebaseline_after_continents: bool,
    /// Enable rigid plate motion (advect plate_id and update velocities)
    pub do_rigid_motion: bool,
    /// Enable collision orogeny (C–C sutures)
    pub do_orogeny: bool,
    /// Enable O–C accretion (arc/forearc growth)
    pub do_accretion: bool,
    /// Enable continental rifting and passive margins
    pub do_rifting: bool,
    /// Enable surface processes (erosion, diffusion, sediment transport/deposition)
    pub do_surface: bool,
    /// Parameter set for surface processes.
    pub surface_params: crate::surface::SurfaceParams,
    /// Cadence: run advection every N steps (>=1). When 1, runs each step.
    pub advection_every: u32,
    /// Cadence for transforms
    pub transforms_every: u32,
    /// Cadence for subduction
    pub subduction_every: u32,
    /// Cadence for flexure
    pub flexure_every: u32,
    /// Cadence for sea-level/isostasy
    pub sea_every: u32,
    /// Gate advection explicitly (combined with cadence)
    pub do_advection: bool,
    /// Gate sea-level explicitly (combined with cadence)
    pub do_sea: bool,
}

/// Result summary for one step.
#[derive(Clone, Copy, Debug)]
pub struct StepStats {
    /// Simulation time after the step (Myr).
    pub t_myr: f64,
    /// Time step size used (Myr).
    pub dt_myr: f64,
    /// Number of divergent boundary edges.
    pub div_count: u32,
    /// Number of convergent boundary edges.
    pub conv_count: u32,
    /// Number of transform boundary edges.
    pub trans_count: u32,
    /// Area-weighted mean of C (0..1).
    pub c_bar: f64,
    /// Flexure residual if applied; otherwise negative.
    pub flex_residual: f32,
}

/// Execute a function for each tile in a `TilePlan` in deterministic order.
/// The closure receives a borrowed `Tile` and a scratch area (caller-owned) for staging.
pub fn for_each_tile<F>(grid: &Grid, plan: &crate::tileplan::TilePlan, mut f: F)
where
    F: FnMut(&crate::tileplan::Tile),
{
    let _ = grid; // grid is provided to mirror planned API and for future use
    for t in &plan.tiles {
        f(t);
    }
}

/// Run the world forward until `t_end_myr` using repeated `step_once` calls.
/// Stepping is chunked into at most `max_steps_per_yield` iterations per call so a
/// caller (e.g., viewer) can interleave UI work between chunks.
pub fn run_to_t(world: &mut World, sp: &StepParams, t_end_myr: f64, max_steps_per_yield: u32) {
    let max_chunk = max_steps_per_yield.max(1);
    while world.clock.t_myr < t_end_myr {
        for _ in 0..max_chunk {
            if world.clock.t_myr >= t_end_myr {
                break;
            }
            let _ = step_once(world, sp);
        }
        // yield to caller
    }
}

/// Execute one evolution step with a minimal CPU pipeline.
///
/// Order:
/// A) age += dt; B) velocities; C) continents advect; D) boundaries classify;
/// E) ridge birth; F) subduction; G) transforms; H) continents uplift;
/// I) flexure; J) isostasy; clock += dt.
pub fn step_once(world: &mut World, sp: &StepParams) -> StepStats {
    // Micro-profiler accumulators (ms)
    let mut ms_boundaries = 0.0f64;
    let mut ms_advection = 0.0f64;
    let mut ms_ridge = 0.0f64;
    let mut ms_subduction = 0.0f64;
    let mut ms_transforms = 0.0f64;
    let mut ms_flexure = 0.0f64;
    let mut ms_sea = 0.0f64;
    // Skip counters for cadence
    let mut sk_advection: u32 = 0;
    let mut sk_transforms: u32 = 0;
    let mut sk_subduction: u32 = 0;
    let mut sk_flexure: u32 = 0;
    let mut sk_sea: u32 = 0;
    let n = world.grid.cells;
    let dt_myr_f32 = sp.dt_myr as f32;

    // B) velocities from plates using current plate ids (or static if disabled)
    let mut vel3 = if sp.do_rigid_motion {
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
        crate::plates::velocity_field_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id)
    } else {
        world.v_en.clone_from(&world.plates.vel_en);
        vec![[0.0f32, 0.0f32, 0.0f32]; world.grid.cells]
    };

    // Diagnostics for velocities
    let mut vmin = f64::INFINITY;
    let mut vmax = 0.0f64;
    let mut vsum = 0.0f64;
    for v in &world.v_en {
        let m = ((v[0] as f64).hypot(v[1] as f64)).abs();
        vmin = vmin.min(m);
        vmax = vmax.max(m);
        vsum += m;
    }
    let vmean = if n > 0 { vsum / (n as f64) } else { 0.0 };

    // A) age simple increment (we will overwrite at ridges below)
    for a in world.age_myr.iter_mut() {
        *a += dt_myr_f32;
    }

    // C) advect continents
    // Determine cadence phase using step_idx (1-based to avoid all-zero edge case)
    let k = world.clock.step_idx.saturating_add(1);
    let do_adv_step = k % sp.advection_every.max(1) as u64 == 0;
    let do_trf_step = sp.do_transforms && (k % (sp.transforms_every.max(1) as u64) == 0);
    let do_sub_step = sp.do_subduction && (k % ((sp.subduction_every.max(1) as u64) * 2) == 0);
    let do_flx_step = sp.do_flexure && (k % (sp.flexure_every.max(1) as u64) == 0);
    let do_sea_step = sp.do_isostasy && (k % sp.sea_every.max(1) as u64 == 0);

    if sp.do_continents && do_adv_step {
        // Seed once if empty (simple deterministic template and uniform thickness)
        if world.c.iter().all(|&c| c == 0.0) && world.th_c_m.iter().all(|&x| x == 0.0) {
            let cp = continent::ContinentParams {
                seed: 1_234_567,
                n_continents: 3,
                mean_radius_km: 2200.0,
                falloff_km: 600.0,
                plateau_uplift_m: 1.0,
                target_land_fraction: None,
            };
            let cf = continent::build_continents(&world.grid, cp);
            // Normalize template into C (0..1)
            world.c.clone_from(&cf.uplift_template_m);
            // Uniform initial thickness 2500 m
            world.th_c_m.fill(2500.0);
        }
        let t0 = Instant::now();
        continent::advect_c_thc(
            &world.grid,
            &world.v_en,
            sp.dt_myr,
            &mut world.c,
            &mut world.th_c_m,
        );
        ms_advection += t0.elapsed().as_secs_f64() * 1000.0;
    } else if sp.do_continents && !do_adv_step {
        sk_advection += 1;
    }

    // C) advect plate_id via semi-Lagrangian nearest-neighbor
    if sp.do_rigid_motion {
        let mut pid_new = world.plates.plate_id.clone();
        crate::sl_advect::advect_plate_id(
            &world.grid,
            &vel3,
            sp.dt_myr,
            &world.plates.plate_id,
            &mut pid_new,
        );
        world.plates.plate_id = pid_new;
        // Refresh velocities after plate-id updates so boundaries and later stages use consistent kinematics
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
        vel3 = crate::plates::velocity_field_m_per_yr(
            &world.grid,
            &world.plates,
            &world.plates.plate_id,
        );
    }

    // D) boundaries classify
    const TAU_OPEN_M_PER_YR: f64 = 0.005;
    let t0b = Instant::now();
    world.boundaries =
        Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, TAU_OPEN_M_PER_YR);
    ms_boundaries += t0b.elapsed().as_secs_f64() * 1000.0;

    // E) ridge birth: reset ages on divergent cells (placeholder: edge-based)
    if sp.do_ridge_birth {
        // Use existing ridge helper on a clone to decide resets
        let t0r = Instant::now();
        let mut ages_tmp = world.age_myr.clone();
        let _ridge_stats = ridge::apply_ridge(
            &world.grid,
            &world.boundaries,
            &mut ages_tmp,
            ridge::RidgeParams { fringe_age_myr: 0.0 },
        );
        for (aw, ar) in world.age_myr.iter_mut().zip(ages_tmp.iter()) {
            if *ar == 0.0 {
                *aw = 0.0;
            }
        }
        ms_ridge += t0r.elapsed().as_secs_f64() * 1000.0;
    }

    // Baseline bathymetry from age (overwrite)
    for i in 0..n {
        let mut d = crate::age::depth_from_age_plate(
            world.age_myr[i] as f64,
            2600.0,
            world.clock.t_myr,
            6000.0,
            1.0e-6,
        ) as f32;
        if !d.is_finite() {
            d = 6000.0;
        }
        world.depth_m[i] = d.clamp(0.0, 6000.0);
    }
    // Prepare to compose tectonic edits as a single delta field added to baseline
    let base_depth = world.depth_m.clone();
    let mut delta_tect: Vec<f32> = vec![0.0; n];

    // F) subduction bands
    let mut last_sub_masks: Option<subduction::SubductionMasks> = None;
    if sp.do_subduction && do_sub_step {
        let t0s = Instant::now();
        let mut sub_tmp = vec![0.0f32; n];
        let sub_res = subduction::apply_subduction(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            &mut sub_tmp,
            subduction::SubductionParams {
                tau_conv_m_per_yr: TAU_OPEN_M_PER_YR,
                trench_half_width_km: 40.0,
                arc_offset_km: 140.0,
                arc_half_width_km: 25.0,
                backarc_width_km: 120.0,
                trench_deepen_m: 1800.0,
                arc_uplift_m: -300.0,
                backarc_uplift_m: -120.0,
                rollback_offset_m: 0.0,
                rollback_rate_km_per_myr: 0.0,
                backarc_extension_mode: false,
                backarc_extension_deepen_m: 400.0,
                continent_c_min: 0.6,
            },
            Some(&world.c),
        );
        for i in 0..n {
            let v = sub_tmp[i];
            if v != 0.0 {
                delta_tect[i] += v - base_depth[i];
            }
        }
        last_sub_masks = Some(sub_res.masks);
        ms_subduction += t0s.elapsed().as_secs_f64() * 1000.0;
    } else if sp.do_subduction && !do_sub_step {
        sk_subduction += 1;
    }

    // F.5) O–C accretion: after subduction, before transforms/orogeny
    if sp.do_accretion {
        // Use masks from last subduction run if available; otherwise recompute quickly with same params
        let sub_masks_owned;
        let sub_masks = if let Some(m) = &last_sub_masks {
            m
        } else {
            let mut tmp = vec![0.0f32; n];
            let sub = subduction::apply_subduction(
                &world.grid,
                &world.boundaries,
                &world.plates.plate_id,
                &world.age_myr,
                &world.v_en,
                &mut tmp,
                subduction::SubductionParams {
                    tau_conv_m_per_yr: TAU_OPEN_M_PER_YR,
                    trench_half_width_km: 40.0,
                    arc_offset_km: 140.0,
                    arc_half_width_km: 25.0,
                    backarc_width_km: 120.0,
                    trench_deepen_m: 1800.0,
                    arc_uplift_m: -300.0,
                    backarc_uplift_m: -120.0,
                    rollback_offset_m: 0.0,
                    rollback_rate_km_per_myr: 0.0,
                    backarc_extension_mode: false,
                    backarc_extension_deepen_m: 400.0,
                    continent_c_min: 0.6,
                },
                Some(&world.c),
            );
            sub_masks_owned = sub.masks;
            &sub_masks_owned
        };
        let p_acc = crate::accretion::AccretionParams {
            k_arc: 0.002,
            gamma_obliquity: 1.0,
            beta_arc: 0.02,
            alpha_arc: 0.001,
            alpha_forearc: 0.0004,
            c_min_continent: 0.6,
            thc_min_m: 0.0,
            thc_max_m: 70_000.0,
            enable_docking: false,
            c_terrane_min: 0.5,
            d_dock_km: 150.0,
            vn_min_m_per_yr: 0.005,
            tau_dock: 0.015,
            couple_flexure: false,
        };
        let mut acc_tmp = vec![0.0f32; n];
        let _ = crate::accretion::apply_oc_accretion(
            &world.grid,
            sub_masks,
            &world.boundaries,
            &vel3,
            &world.plates.kind,
            &mut world.c,
            &mut world.th_c_m,
            &mut acc_tmp,
            &world.area_m2,
            &p_acc,
            sp.dt_myr,
        );
        for i in 0..n {
            delta_tect[i] += acc_tmp[i];
        }
    }

    // F.8) Continental rifting before transforms/orogeny/flexure
    if sp.do_rifting {
        let p_rift = crate::rifting::RiftingParams {
            c_rift_min: 0.6,
            v_open_min_m_per_yr: 0.001,
            w_core_km: 60.0,
            w_taper_km: 250.0,
            k_thin: 0.10,
            alpha_subs: 0.8,
            ocean_thresh: 0.15,
            k_c_oceanize: 0.03,
            reset_age_on_core: true,
            enable_shoulder: false,
            w_bulge_km: 120.0,
            beta_shoulder: 0.2,
            couple_flexure: false,
            thc_min_m: 20_000.0,
            thc_max_m: 70_000.0,
        };
        let mut rift_tmp = vec![0.0f32; n];
        let _ = crate::rifting::apply_rifting(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.plates.kind,
            &mut world.c,
            &mut world.th_c_m,
            &mut world.age_myr,
            &mut rift_tmp,
            &world.area_m2,
            &p_rift,
            sp.dt_myr,
        );
        for i in 0..n {
            delta_tect[i] += rift_tmp[i];
        }
    }

    // G) transforms
    if sp.do_transforms && do_trf_step {
        let t0t = Instant::now();
        let mut tr_tmp = vec![0.0f32; n];
        let _ = transforms::apply_transforms(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.v_en,
            &mut tr_tmp,
            transforms::TransformParams {
                tau_open_m_per_yr: TAU_OPEN_M_PER_YR,
                min_tangential_m_per_yr: 0.003,
                max_normal_m_per_yr: 0.010,
                basin_half_width_km: 50.0,
                ridge_like_uplift_m: -200.0,
                basin_deepen_m: 300.0,
            },
            sp.dt_myr,
        );
        for i in 0..n {
            delta_tect[i] += tr_tmp[i];
        }
        ms_transforms += t0t.elapsed().as_secs_f64() * 1000.0;
    } else if sp.do_transforms && !do_trf_step {
        sk_transforms += 1;
    }

    // G.5) Orogeny (C–C): after transforms, before continents uplift/flexure
    if sp.do_orogeny {
        let p_orog = crate::orogeny::OrogenyParams {
            c_min: 0.6,
            w_core_km: 120.0,
            w_taper_km: 220.0,
            k_thick: 0.08,
            beta_uplift: 0.8,
            gamma_obliquity: 1.0,
            couple_flexure: false,
        };
        let mut orog_tmp = vec![0.0f32; n];
        let _stats = crate::orogeny::apply_cc_orogeny(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.c,
            &world.area_m2,
            &mut world.th_c_m,
            &mut orog_tmp,
            &p_orog,
            sp.dt_myr,
            &world.plates.kind,
        );
        for i in 0..n {
            delta_tect[i] += orog_tmp[i];
        }
    }

    // H) continents uplift
    if sp.do_continents {
        // NOTE: This applies uplift cumulatively each step as `depth += -(C * th_c)`, independent
        // of `dt`. If called every step, depths can drift unbounded. Consider computing a
        // diagnostic uplift field and composing it (once) with the age baseline instead of adding
        // per-step, or scale by `dt` and clamp to a target uplift amplitude tied to mass balance.
        let mut cont_tmp = vec![0.0f32; n];
        continent::apply_uplift_from_c_thc(&mut cont_tmp, &world.c, &world.th_c_m);
        for i in 0..n {
            delta_tect[i] += cont_tmp[i];
        }
        // Consider this a continents epoch change only if the uplift changed depths appreciably
        // (we simply bump every time continents are applied in this MVP)
        world.epoch_continents = world.epoch_continents.wrapping_add(1);
    }

    // Finalize tectonic composition: depth = baseline + accumulated deltas
    for i in 0..n {
        let v = base_depth[i] + delta_tect[i];
        world.depth_m[i] = v.clamp(-8000.0, 8000.0);
    }

    // H.5) Surface processes after tectonics, before flexure/isostasy
    world.last_surface_stats = None;
    // If flexure should be applied after surface coupling, skip the later generic flexure stage
    let mut run_flexure_stage_i = sp.do_flexure;
    if sp.do_surface {
        let stats = crate::surface::apply_surface_processes(
            &world.grid,
            &world.c,
            &mut world.depth_m,
            &mut world.sediment_m,
            &world.area_m2,
            &sp.surface_params,
            sp.dt_myr,
        );
        // Log once per step when enabled
        let res_pct = if stats.eroded_m3 > 0.0 { stats.residual_m3 / stats.eroded_m3 } else { 0.0 };
        println!(
            "[surface] eroded={:.2e} m³ deposited={:.2e} m³ residual={:+.3}% max_ero={:.2} m max_dep={:.2} m",
            stats.eroded_m3, stats.deposited_m3, res_pct * 100.0, stats.max_erosion_m, stats.max_deposition_m
        );
        world.last_surface_stats = Some(stats);
        if sp.surface_params.couple_flexure {
            // Run a quick Winkler-like flexure response to approximate coupling
            world.last_flex_residual = -1.0;
            let tf0 = Instant::now();
            let lp = flexure_loads::LoadParams {
                rho_w: 1030.0,
                rho_c: 2900.0,
                g: 9.81,
                sea_level_m: world.sea.eta_m,
            };
            let f_load = flexure_loads::assemble_load_from_depth(&world.grid, &world.depth_m, &lp);
            let k = 3.0e8f32;
            let mut w = vec![0.0f32; n];
            for i in 0..n {
                w[i] = f_load[i] / k;
            }
            let mut r0: f64 = 0.0;
            let mut r1: f64 = 0.0;
            for i in 0..n {
                let fi = f_load[i] as f64;
                r0 += fi * fi;
                let ri = fi - (k as f64) * (w[i] as f64);
                r1 += ri * ri;
            }
            world.last_flex_residual = ((r1.sqrt()) / (r0.sqrt().max(1.0))) as f32;
            for (d, wi) in world.depth_m.iter_mut().zip(w.iter()) {
                *d = (*d + *wi).clamp(-8000.0, 8000.0);
            }
            // Skip generic flexure stage below if we already coupled here
            run_flexure_stage_i = false;
            ms_flexure += tf0.elapsed().as_secs_f64() * 1000.0;
        }
    }

    // I) flexure (Winkler placeholder)
    world.last_flex_residual = -1.0;
    if run_flexure_stage_i && do_flx_step {
        let tf1 = Instant::now();
        let lp = flexure_loads::LoadParams {
            rho_w: 1030.0,
            rho_c: 2900.0,
            g: 9.81,
            sea_level_m: world.sea.eta_m,
        };
        let f_load = flexure_loads::assemble_load_from_depth(&world.grid, &world.depth_m, &lp);
        let k = 3.0e8f32;
        let mut w = vec![0.0f32; n];
        for i in 0..n {
            w[i] = f_load[i] / k;
        }
        let mut r0: f64 = 0.0;
        let mut r1: f64 = 0.0;
        for i in 0..n {
            let fi = f_load[i] as f64;
            r0 += fi * fi;
            let ri = fi - (k as f64) * (w[i] as f64);
            r1 += ri * ri;
        }
        r0 = r0.sqrt();
        r1 = r1.sqrt();
        world.last_flex_residual = (r1 / r0.max(1.0)) as f32;
        for (d, wi) in world.depth_m.iter_mut().zip(w.iter()) {
            *d = (*d + *wi).clamp(-8000.0, 8000.0);
        }
        ms_flexure += tf1.elapsed().as_secs_f64() * 1000.0;
    } else if run_flexure_stage_i && !do_flx_step {
        sk_flexure += 1;
    }

    // J) isostasy / sea-level: maintain reference ocean volume if available, plus slow eustasy
    if sp.do_isostasy && do_sea_step {
        let ts0 = Instant::now();
        if world.sea_level_ref.is_none() {
            world.sea_level_ref =
                Some(isostasy::compute_ref(&world.depth_m, &world.area_m2, world.sea.eta_m));
        }
        if sp.auto_rebaseline_after_continents
            && world.epoch_continents != world.last_rebaseline_epoch
        {
            let area_clone: Vec<f32> = world.area_m2.clone();
            let _ = isostasy::rebaseline(world, &area_clone);
            world.last_rebaseline_epoch = world.epoch_continents;
        }
        if let Some(r) = world.sea_level_ref {
            // Isostatic offset for target volume
            let l_iso = isostasy::solve_offset_for_volume(
                &world.depth_m,
                &world.area_m2,
                r.volume_m3,
                1e6,
                64,
            );
            // Eustasy extra offset (policy: constant 0 in MVP)
            let policy = crate::sea_level::EustasyPolicy::Constant { eta_m: 0.0 };
            let eta = crate::sea_level::update_eustasy_eta(
                world.clock.t_myr,
                sp.dt_myr,
                &policy,
                &world.depth_m,
                &world.area_m2,
                l_iso,
            );
            let l_total = l_iso + eta;
            // Unified sea-level policy: store in `world.sea.eta_m` and do not mutate `depth_m`.
            world.sea.eta_m = l_total as f32;
            let ocean_frac = {
                let mut ocean_area = 0.0f64;
                let mut total_area = 0.0f64;
                for i in 0..n {
                    // Evaluate ocean with eta applied: depth + eta > 0
                    if (world.depth_m[i] as f64 + l_total) > 0.0 {
                        ocean_area += world.area_m2[i] as f64;
                    }
                    total_area += world.area_m2[i] as f64;
                }
                if total_area > 0.0 {
                    ocean_area / total_area
                } else {
                    0.0
                }
            };
            println!(
                "[sea] L_iso={:+.1} m  eta={:+.1} m  policy=constant  ocean={:.1}%",
                l_iso,
                eta,
                ocean_frac * 100.0
            );
            if l_iso.abs() < 1e-3 {
                println!("[isostasy] offset≈0 after rebaseline (|Δ|={:.4} m)", l_iso.abs());
            }
        }
        ms_sea += ts0.elapsed().as_secs_f64() * 1000.0;
    } else if sp.do_isostasy && !do_sea_step {
        sk_sea += 1;
    }

    // Update clock
    world.clock.t_myr += sp.dt_myr;
    world.clock.step_idx = world.clock.step_idx.saturating_add(1);

    // Stats for log line
    let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
    let mut weighted_c: f64 = 0.0;
    for i in 0..n {
        weighted_c += (world.c[i] as f64) * (world.area_m2[i] as f64);
    }
    let c_bar = if total_area > 0.0 { weighted_c / total_area } else { 0.0 };

    // Perf log
    let step_ms = ms_boundaries
        + ms_advection
        + ms_ridge
        + ms_subduction
        + ms_transforms
        + ms_flexure
        + ms_sea;
    println!(
        "[perf] step={:.2} ms | boundaries={:.2} | sea={:.2}/skip:{} | advection={:.2}/skip:{} | transforms={:.2}/skip:{} | subduction={:.2}/skip:{} | flexure={:.2}/skip:{}",
        step_ms,
        ms_boundaries,
        ms_sea,
        sk_sea,
        ms_advection,
        sk_advection,
        ms_transforms,
        sk_transforms,
        ms_subduction,
        sk_subduction,
        ms_flexure,
        sk_flexure
    );

    // Advect diagnostics (rough backtrace distance)
    let backtrace_km = (vmax * (sp.dt_myr * 1.0e6)) / 1000.0;
    println!(
        "[advect] |V| m/yr min/mean/max = {:.6}/{:.6}/{:.6}  backtrace max ≈ {:.1} km dt={:.1} Myr",
        if vmin.is_finite() { vmin } else { 0.0 },
        vmean,
        vmax,
        backtrace_km,
        sp.dt_myr
    );
    if sp.do_continents {
        let mut thc_sum = 0.0f64;
        for &x in &world.th_c_m {
            thc_sum += x as f64;
        }
        let thc_mean_km = if n > 0 { (thc_sum / (n as f64)) / 1000.0 } else { 0.0 };
        println!("[continents] C̄={:.1}% th_c mean={:.2} km", c_bar * 100.0, thc_mean_km);
    }

    StepStats {
        t_myr: world.clock.t_myr,
        dt_myr: sp.dt_myr,
        div_count: world.boundaries.stats.divergent,
        conv_count: world.boundaries.stats.convergent,
        trans_count: world.boundaries.stats.transform,
        c_bar,
        flex_residual: world.last_flex_residual,
    }
}
