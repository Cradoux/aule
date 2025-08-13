//! World state container and constructors (T-400).

use crate::{
    age::depth_from_age, boundaries::Boundaries, continent, flexure_loads, grid::Grid, isostasy,
    plates::Plates, ridge, subduction, transforms,
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
    /// Bathymetry depth in meters (+ down), length = grid.cells.
    pub depth_m: Vec<f32>,
    /// Continental fraction (0..1). Seeded when continents are enabled.
    pub c: Vec<f32>,
    /// Continental crust thickness in meters (≥0). Seeded alongside `C`.
    pub th_c_m: Vec<f32>,
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
    /// Last surface-processes stats (if applied this step).
    pub last_surface_stats: Option<crate::surface::SurfaceStats>,
    /// Continents change epoch counter (bump when C/th_c change applied).
    pub epoch_continents: u64,
    /// Last epoch when sea-level was re-baselined (to debounce auto mode).
    pub last_rebaseline_epoch: u64,
}

/// Reference sea-level bookkeeping
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeaLevelRef {
    /// Target ocean volume at reference state (m^3)
    pub volume_m3: f64,
    /// Ocean area used at reference state (m^2)
    pub ocean_area_m2: f64,
}

impl World {
    /// Construct a world with deterministic plates and zeroed ages/depths.
    pub fn new(f: u32, num_plates: u32, seed: u64) -> Self {
        let grid = Grid::new(f);
        let plates = Plates::new(&grid, num_plates, seed);
        let v_en = plates.vel_en.clone();
        let boundaries =
            crate::boundaries::Boundaries::classify(&grid, &plates.plate_id, &v_en, 0.005);
        let age_myr = vec![0.0f32; grid.cells];
        let depth_m = vec![0.0f32; grid.cells];
        let c = vec![0.0f32; grid.cells];
        let th_c_m = vec![0.0f32; grid.cells];
        let sediment_m = vec![0.0f32; grid.cells];
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
            depth_m,
            c,
            th_c_m,
            v_en,
            clock,
            sea_level_ref: None,
            area_m2,
            last_flex_residual: 0.0,
            sediment_m,
            last_surface_stats: None,
            epoch_continents: 0,
            last_rebaseline_epoch: 0,
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
    let n = world.grid.cells;
    let dt_myr_f32 = sp.dt_myr as f32;

    // B) velocities from plates using current plate ids (or static if disabled)
    let vel3 = if sp.do_rigid_motion {
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
    if sp.do_continents {
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
        continent::advect_c_thc(
            &world.grid,
            &world.v_en,
            sp.dt_myr,
            &mut world.c,
            &mut world.th_c_m,
        );
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
    }

    // D) boundaries classify
    const TAU_OPEN_M_PER_YR: f64 = 0.005;
    world.boundaries =
        Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, TAU_OPEN_M_PER_YR);

    // E) ridge birth: reset ages on divergent cells (placeholder: edge-based)
    if sp.do_ridge_birth {
        // Use existing ridge helper on a clone to decide resets
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
    }

    // Baseline bathymetry from age (overwrite)
    for i in 0..n {
        let mut d = depth_from_age(world.age_myr[i] as f64, 2600.0, 350.0, 0.0) as f32;
        if !d.is_finite() {
            d = 6000.0;
        }
        world.depth_m[i] = d.clamp(0.0, 6000.0);
    }

    // F) subduction bands
    if sp.do_subduction {
        let _sub = subduction::apply_subduction(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            &mut world.depth_m,
            subduction::SubductionParams {
                tau_conv_m_per_yr: TAU_OPEN_M_PER_YR,
                trench_half_width_km: 50.0,
                arc_offset_km: 150.0,
                arc_half_width_km: 30.0,
                backarc_width_km: 150.0,
                trench_deepen_m: 3000.0,
                arc_uplift_m: -500.0,
                backarc_uplift_m: -200.0,
                rollback_offset_m: 0.0,
                rollback_rate_km_per_myr: 0.0,
                backarc_extension_mode: false,
                backarc_extension_deepen_m: 600.0,
            },
        );
    }

    // F.5) O–C accretion: after subduction, before transforms/orogeny
    if sp.do_accretion {
        // Use masks from last subduction run if available by recomputing quickly with same params
        let sub = subduction::apply_subduction(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            &mut world.depth_m,
            subduction::SubductionParams {
                tau_conv_m_per_yr: TAU_OPEN_M_PER_YR,
                trench_half_width_km: 50.0,
                arc_offset_km: 150.0,
                arc_half_width_km: 30.0,
                backarc_width_km: 150.0,
                trench_deepen_m: 3000.0,
                arc_uplift_m: -500.0,
                backarc_uplift_m: -200.0,
                rollback_offset_m: 0.0,
                rollback_rate_km_per_myr: 0.0,
                backarc_extension_mode: false,
                backarc_extension_deepen_m: 600.0,
            },
        );
        let p_acc = crate::accretion::AccretionParams {
            k_arc: 0.05,
            gamma_obliquity: 1.0,
            beta_arc: 1.0,
            alpha_arc: 0.02,
            alpha_forearc: 0.01,
            c_min_continent: 0.6,
            thc_min_m: 20_000.0,
            thc_max_m: 70_000.0,
            enable_docking: false,
            c_terrane_min: 0.5,
            d_dock_km: 150.0,
            vn_min_m_per_yr: 0.005,
            tau_dock: 0.02,
            couple_flexure: false,
        };
        let _ = crate::accretion::apply_oc_accretion(
            &world.grid,
            &sub.masks,
            &world.boundaries,
            &vel3,
            &mut world.c,
            &mut world.th_c_m,
            &mut world.depth_m,
            &world.area_m2,
            &p_acc,
            sp.dt_myr,
        );
    }

    // F.8) Continental rifting before transforms/orogeny/flexure
    if sp.do_rifting {
        let p_rift = crate::rifting::RiftingParams {
            c_rift_min: 0.6,
            v_open_min_m_per_yr: 0.001,
            w_core_km: 60.0,
            w_taper_km: 250.0,
            k_thin: 0.15,
            alpha_subs: 1.0,
            ocean_thresh: 0.15,
            k_c_oceanize: 0.05,
            reset_age_on_core: true,
            enable_shoulder: false,
            w_bulge_km: 120.0,
            beta_shoulder: 0.3,
            couple_flexure: false,
            thc_min_m: 20_000.0,
            thc_max_m: 70_000.0,
        };
        let _ = crate::rifting::apply_rifting(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &mut world.c,
            &mut world.th_c_m,
            &mut world.age_myr,
            &mut world.depth_m,
            &world.area_m2,
            &p_rift,
            sp.dt_myr,
        );
    }

    // G) transforms
    if sp.do_transforms {
        let _ = transforms::apply_transforms(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.v_en,
            &mut world.depth_m,
            transforms::TransformParams {
                tau_open_m_per_yr: TAU_OPEN_M_PER_YR,
                min_tangential_m_per_yr: 0.002,
                basin_half_width_km: 60.0,
                ridge_like_uplift_m: -300.0,
                basin_deepen_m: 400.0,
            },
        );
    }

    // G.5) Orogeny (C–C): after transforms, before continents uplift/flexure
    if sp.do_orogeny {
        let p_orog = crate::orogeny::OrogenyParams {
            c_min: 0.6,
            w_core_km: 150.0,
            w_taper_km: 250.0,
            k_thick: 0.10,
            beta_uplift: 1.0,
            gamma_obliquity: 1.0,
            couple_flexure: false,
        };
        let _stats = crate::orogeny::apply_cc_orogeny(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.c,
            &world.area_m2,
            &mut world.th_c_m,
            &mut world.depth_m,
            &p_orog,
            sp.dt_myr,
        );
    }

    // H) continents uplift
    if sp.do_continents {
        continent::apply_uplift_from_c_thc(&mut world.depth_m, &world.c, &world.th_c_m);
        // Consider this a continents epoch change only if the uplift changed depths appreciably
        // (we simply bump every time continents are applied in this MVP)
        world.epoch_continents = world.epoch_continents.wrapping_add(1);
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
            let lp = flexure_loads::LoadParams {
                rho_w: 1030.0,
                rho_c: 2900.0,
                g: 9.81,
                sea_level_m: 0.0,
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
        }
    }

    // I) flexure (Winkler placeholder)
    world.last_flex_residual = -1.0;
    if run_flexure_stage_i {
        let lp =
            flexure_loads::LoadParams { rho_w: 1030.0, rho_c: 2900.0, g: 9.81, sea_level_m: 0.0 };
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
    }

    // J) isostasy / sea-level: maintain reference ocean volume if available, plus slow eustasy
    if sp.do_isostasy {
        if world.sea_level_ref.is_none() {
            world.sea_level_ref = Some(isostasy::compute_ref(&world.depth_m, &world.area_m2));
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
            isostasy::apply_sea_level_offset(&mut world.depth_m, l_total);
            let ocean_frac = {
                let mut ocean_area = 0.0f64;
                let mut total_area = 0.0f64;
                for i in 0..n {
                    if world.depth_m[i] > 0.0 {
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
