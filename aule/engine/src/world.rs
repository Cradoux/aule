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
    boundaries::Boundaries, grid::Grid, isostasy, plates::Plates,
};

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
        let pc = crate::PhysConsts::default();
        let r_earth_m = pc.r_earth_m;
        let scale = 4.0 * std::f64::consts::PI * r_earth_m * r_earth_m;
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
        // CRITICAL: Ensure all plates start with zero angular velocity to prevent persistence
        for omega in &mut self.plates.omega_rad_yr {
            *omega = 0.0;
        }
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