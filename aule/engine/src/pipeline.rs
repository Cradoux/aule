//! Centralized physics pipeline (T-632): run full step sequence and expose a single surface elevation.
//!
//! This viewer-facing pipeline contrasts with `world::step_once` in one key way:
//! - Here we compute and output `eta` (sea level) but do not bake the offset back into `world.depth_m`.
//!   Callers must render elevation as `z = -depth - eta`. This reduces coastline flicker and keeps
//!   `depth_m` as a separable tectonic + age baseline.
//! - Subduction/transform/rifting/orogeny still edit `depth_m` additively on top of age-derived bathy.
//!   If these are applied at high cadence without normalization, long-run drift of `depth_m` can occur.

use crate::{boundaries, flexure_loads, ridge, subduction, transforms, world::World};

/// Pipeline configuration used by Simple and Advanced modes.
#[derive(Clone, Copy, Debug)]
pub struct PipelineCfg {
    /// Time step in Myr.
    pub dt_myr: f32,
    /// Number of steps to run per frame (viewer may loop externally as well).
    pub steps_per_frame: u32,
    /// Enable flexure solve stage.
    pub enable_flexure: bool,
    /// Enable erosion/diffusion stage.
    pub enable_erosion: bool,
    /// Target land fraction (0..1) used by sea-level solve.
    pub target_land_frac: f32,
    /// If true, keep eta fixed (skip sea-level regulation this tick).
    pub freeze_eta: bool,
    /// If true, append a one-line mass budget log each step.
    pub log_mass_budget: bool,
    /// Enable subduction band edits (disable to isolate physics noise)
    pub enable_subduction: bool,
    /// Enable rigid plate motion (advect `plate_id`, refresh velocities)
    pub enable_rigid_motion: bool,
    /// Cadence controls (>=1); a stage runs when `(step_idx+1) % cadence == 0`.
    /// Transformer bands cadence in steps.
    pub cadence_trf_every: u32,
    /// Subduction bands cadence in steps.
    pub cadence_sub_every: u32,
    /// Flexure cadence in steps.
    pub cadence_flx_every: u32,
    /// Sea-level (eta) solve cadence in steps.
    pub cadence_sea_every: u32,
    /// Surface processes cadence in steps.
    pub cadence_surf_every: u32,
    /// Sub-steps for transforms per cadence (narrow operator stability)
    pub substeps_transforms: u32,
    /// Sub-steps for subduction band deltas per cadence
    pub substeps_subduction: u32,
    /// If true, attempt GPU flexure (experimental). Fallback to Winkler if unavailable.
    pub use_gpu_flexure: bool,
    /// GPU flexure: number of multigrid levels.
    pub gpu_flex_levels: u32,
    /// GPU flexure: V-cycles per cadence.
    pub gpu_flex_cycles: u32,
    /// GPU flexure: weighted-Jacobi omega.
    pub gpu_wj_omega: f32,
    /// Subtract mean load before solving (stability, removes rigid-body mode).
    pub subtract_mean_load: bool,
    // Subduction tunables (viewer-controlled)
    /// Convergence threshold in m/yr used for band detection.
    pub sub_tau_conv_m_per_yr: f32,
    /// Trench band half-width in km on subducting side.
    pub sub_trench_half_width_km: f32,
    /// Arc band peak offset from trench hinge in km (overriding side).
    pub sub_arc_offset_km: f32,
    /// Arc band half-width in km.
    pub sub_arc_half_width_km: f32,
    /// Back-arc band width in km (behind arc).
    pub sub_backarc_width_km: f32,
    /// Trench deepening magnitude in meters (positive deepens).
    pub sub_trench_deepen_m: f32,
    /// Arc uplift magnitude in meters (negative uplifts/shallows).
    pub sub_arc_uplift_m: f32,
    /// Back-arc uplift magnitude in meters (negative uplifts/shallows).
    pub sub_backarc_uplift_m: f32,
    /// Absolute rollback offset applied to overriding distances in meters.
    pub sub_rollback_offset_m: f32,
    /// Rollback rate in km/Myr for time-progression (applied by viewer logic).
    pub sub_rollback_rate_km_per_myr: f32,
    /// If true, back-arc is deepened (extension mode) instead of uplifted.
    pub sub_backarc_extension_mode: bool,
    /// Back-arc deepening magnitude in meters when in extension mode.
    pub sub_backarc_extension_deepen_m: f32,
    /// Continental fraction threshold [0,1] for gating trench deepening.
    pub sub_continent_c_min: f32,

    // Surface processes parameters (viewer-controlled)
    /// Stream-power coefficient for fluvial incision (yr⁻¹ m^(1-m) s^m units folded)
    pub surf_k_stream: f32,
    /// Stream-power m exponent
    pub surf_m_exp: f32,
    /// Stream-power n exponent
    pub surf_n_exp: f32,
    /// Hillslope diffusion coefficient κ (m²/yr)
    pub surf_k_diff: f32,
    /// Sediment transport coefficient (nondimensional scaling)
    pub surf_k_tr: f32,
    /// Transport-limited p exponent
    pub surf_p_exp: f32,
    /// Transport-limited q exponent
    pub surf_q_exp: f32,
    /// Sediment density (kg/m³)
    pub surf_rho_sed: f32,
    /// Minimum slope threshold
    pub surf_min_slope: f32,
    /// Subcycles per step for numerical stability (>=1)
    pub surf_subcycles: u32,
    /// If true, couple surface mass change back into flexure load (placeholder)
    pub surf_couple_flexure: bool,

    // Plate lifecycle (scaffolding)
    /// If >0, every N steps attempt to spawn a microplate at a random convergent cell.
    pub cadence_spawn_plate_every: u32,
    /// If >0, every N steps attempt to retire the smallest-degree plate into a neighbor.
    pub cadence_retire_plate_every: u32,

    // Force-balance Euler update
    /// If >0, every N steps apply a small force-balance update to plate omegas.
    pub cadence_force_balance_every: u32,
    /// Force-balance gain (unitless scale to map boundary drive to dω).
    pub fb_gain: f32,
    /// Force-balance damping per Myr.
    pub fb_damp_per_myr: f32,
    /// Weight for convergent edges.
    pub fb_k_conv: f32,
    /// Weight for divergent edges.
    pub fb_k_div: f32,
    /// Weight for transform edges.
    pub fb_k_trans: f32,
    /// Max per-step Δω clamp.
    pub fb_max_domega: f32,
    /// Max absolute |ω| clamp.
    pub fb_max_omega: f32,
}

/// Borrowed views of surface fields required by the pipeline.
pub struct SurfaceFields<'a> {
    /// Solved surface elevation (meters), vertex-major, same order as GPU uploads.
    pub elev_m: &'a mut [f32],
    /// Sea level eta (meters) relative to geoid; output of sea-level solve.
    pub eta_m: &'a mut f32,
}

/// Execute the full physics step sequence once.
///
/// Order:
/// 1) boundaries
/// 2) advection (kinematics)
/// 3) ridges (oceanic crust production)
/// 4) subduction (convergence thinning/consumption)
/// 5) flexure (optional)
/// 6) erosion/diffusion (optional)
/// 7) sea-level bisection to target land fraction
/// 8) sanitize elevation
pub fn step_full(world: &mut World, surf: SurfaceFields, cfg: PipelineCfg) {
    let dt = cfg.dt_myr;
    let k = world.clock.step_idx.saturating_add(1);
    let do_trf = k % (cfg.cadence_trf_every.max(1) as u64) == 0;
    let do_sub = k % (cfg.cadence_sub_every.max(1) as u64) == 0;
    let do_flx = k % (cfg.cadence_flx_every.max(1) as u64) == 0;
    let do_sea = k % (cfg.cadence_sea_every.max(1) as u64) == 0;
    let do_spawn =
        cfg.cadence_spawn_plate_every > 0 && k % (cfg.cadence_spawn_plate_every as u64) == 0;
    let do_retire =
        cfg.cadence_retire_plate_every > 0 && k % (cfg.cadence_retire_plate_every as u64) == 0;
    let do_fb =
        cfg.cadence_force_balance_every > 0 && k % (cfg.cadence_force_balance_every as u64) == 0;

    // Kinematics:
    // - Refresh per-cell velocities from plates
    // - Optionally advect plate_id via semi-Lagrangian nearest-neighbour when enabled
    // - Classify boundaries using the refreshed fields
    // Accumulate 3D velocities for processes that need them
    let mut vel3: Vec<[f32; 3]> = vec![[0.0; 3]; world.grid.cells];
    if cfg.enable_rigid_motion {
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
        vel3 = crate::plates::velocity_field_m_per_yr(
            &world.grid,
            &world.plates,
            &world.plates.plate_id,
        );
        // Subcycle semi-Lagrangian advection if displacement per step exceeds a fraction of mean cell size
        let dt_myr = cfg.dt_myr as f64;
        let mut vmax_m_per_yr: f64 = 0.0;
        for v in &vel3 {
            let m = ((v[0] as f64).hypot(v[1] as f64)).hypot(v[2] as f64);
            if m > vmax_m_per_yr {
                vmax_m_per_yr = m;
            }
        }
        // Estimate mean linear cell scale from area
        let mut area_mean: f64 = 0.0;
        for a in &world.area_m2 {
            area_mean += *a as f64;
        }
        let n_cells = world.grid.cells.max(1) as f64;
        area_mean /= n_cells;
        let cell_scale_m = area_mean.max(1.0).sqrt();
        let disp = vmax_m_per_yr * (dt_myr * 1.0e6);
        // Aggressive subcycling so each substep moves < 0.1 cell-scale
        let substeps = ((disp / (0.1 * cell_scale_m)).ceil() as u32).clamp(1, 32);
        let sub_dt = dt_myr / (substeps as f64);
        let mut pid_curr = world.plates.plate_id.clone();
        let mut pid_next = pid_curr.clone();
        for _ in 0..substeps {
            crate::sl_advect::advect_plate_id(&world.grid, &vel3, sub_dt, &pid_curr, &mut pid_next);
            std::mem::swap(&mut pid_curr, &mut pid_next);
        }
        world.plates.plate_id = pid_curr;
        // Post-fix tiny holes: if a cell disagrees with strict majority of its 1-ring neighbors, adopt the majority label
        let mut pid_fixed = world.plates.plate_id.clone();
        for (i, pid_out) in pid_fixed.iter_mut().enumerate().take(world.grid.cells) {
            let mut counts: std::collections::HashMap<u16, u32> = std::collections::HashMap::new();
            for &n in &world.grid.n1[i] {
                let id = world.plates.plate_id[n as usize];
                *counts.entry(id).or_insert(0) += 1;
            }
            if counts.is_empty() {
                continue;
            }
            if let Some((mode_id, mode_count)) =
                counts.iter().max_by_key(|e| e.1).map(|(k, v)| (*k, *v))
            {
                let total = counts.values().sum::<u32>();
                if mode_count > total / 2 && world.plates.plate_id[i] != mode_id {
                    *pid_out = mode_id;
                }
            }
        }
        world.plates.plate_id = pid_fixed;
    } else {
        world.v_en.clone_from(&world.plates.vel_en);
        vel3.fill([0.0, 0.0, 0.0]);
    }

    // Boundaries classification
    world.boundaries =
        boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
    // Plate-id diagnostics and healing
    {
        let mut invalid = 0usize;
        for &pid in &world.plates.plate_id {
            if pid == crate::plates::INVALID_PLATE_ID {
                invalid += 1;
            }
        }
        if invalid > 0 {
            println!("[debug] invalid_plate_cells={}", invalid);
            crate::plates::heal_plate_ids(&world.grid, &mut world.plates.plate_id);
            // Refresh velocities after healing to keep kinematics consistent
            world.v_en = crate::plates::velocity_en_m_per_yr(
                &world.grid,
                &world.plates,
                &world.plates.plate_id,
            );
        }
    }
    // Plate birth/death scaffolding (diagnostic only; disabled by default via cadence=0)
    if do_spawn {
        // Pick a convergent boundary cell if available; fallback to center cell
        let mut pick: Option<usize> = None;
        for &(u, _v, cls) in &world.boundaries.edges {
            if cls == 2 {
                pick = Some(u as usize);
                break;
            }
        }
        let cell = pick.unwrap_or(world.grid.cells / 2);
        let new_id = world.plates.spawn_plate_at(&world.grid, cell, 0xB17F_1A7E, 2);
        println!("[plates] spawn id={} at cell={} (ring=2)", new_id, cell);
        // Refresh velocities used downstream
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
    }
    if do_retire {
        // Retire the highest-degree plate into one of its neighbors (if any), else no-op
        let adj = crate::plate_network::build_plate_adjacency(&world.boundaries, &world.plates);
        let mut best_pid: Option<u16> = None;
        let mut best_deg: i32 = -1;
        for (&pid, nbrs) in &adj {
            let deg = nbrs.len() as i32;
            if deg > best_deg {
                best_deg = deg;
                best_pid = Some(pid);
            }
        }
        if let Some(pid) = best_pid {
            if let Some(nbrs) = adj.get(&pid) {
                if let Some(&to) = nbrs.iter().next() {
                    world.plates.retire_plate_soft(&world.grid, pid, to);
                    println!("[plates] retire id={} → {} (soft)", pid, to);
                    world.v_en = crate::plates::velocity_en_m_per_yr(
                        &world.grid,
                        &world.plates,
                        &world.plates.plate_id,
                    );
                }
            }
        }
    }
    // Plate network diagnostics (lightweight)
    if k % 32 == 0 {
        let net =
            crate::plate_network::compute_stats(&world.grid, &world.plates, &world.boundaries);
        println!("[network] {}", crate::plate_network::summarize(&net));
    }
    if k % 64 == 0 {
        let bh = crate::plate_network::compute_boundary_health(&world.boundaries);
        println!("[network] {}", crate::plate_network::summarize_boundary(&bh));
    }

    // Optional force-balance Euler poles update
    if do_fb {
        let p = crate::force_balance::FbParams {
            gain: cfg.fb_gain,
            damp_per_myr: cfg.fb_damp_per_myr,
            k_conv: cfg.fb_k_conv,
            k_div: cfg.fb_k_div,
            k_trans: cfg.fb_k_trans,
            max_domega: cfg.fb_max_domega,
            max_omega: cfg.fb_max_omega,
        };
        let pid_copy = world.plates.plate_id.clone();
        crate::force_balance::apply_force_balance(
            &world.grid,
            &world.boundaries,
            &pid_copy,
            &mut world.plates,
            &world.area_m2,
            cfg.dt_myr as f64,
            p,
        );
        // Refresh velocities after adjusting omegas
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
    }

    // Age–Depth residual and hypsometry diagnostics (lightweight)
    if k % 16 == 0 {
        // Age→depth predictor using plate-cooling blend (same params as world gen)
        let d_pred: Vec<f32> = world
            .age_myr
            .iter()
            .map(|&age| {
                crate::age::depth_from_age_blend(age as f64, 2600.0, 350.0, 6300.0, 60.0) as f32
            })
            .collect();
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        let mut wsum = 0.0f64;
        let mut absmax = 0.0f32;
        for (i, _) in d_pred.iter().enumerate().take(world.grid.cells) {
            let w = world.area_m2[i] as f64;
            let r = world.depth_m[i] - d_pred[i];
            sum += (r as f64) * w;
            sumsq += (r as f64) * (r as f64) * w;
            wsum += w;
            if r.abs() > absmax {
                absmax = r.abs();
            }
        }
        if wsum > 0.0 {
            let mean = sum / wsum;
            let var = (sumsq / wsum) - mean * mean;
            let std = var.max(0.0).sqrt();
            println!(
                "[diag] age-depth residual: mean={:.1} m std={:.1} m max|r|={:.0} m",
                mean, std, absmax
            );
        }
        // Hypsometry snapshot (ocean thickness positive)
        let mut land_area = 0.0f64;
        let mut ocean_area = 0.0f64;
        for i in 0..world.grid.cells {
            let a = world.area_m2[i] as f64;
            let z = (*surf.eta_m as f64) - (world.depth_m[i] as f64);
            if z > 0.0 {
                land_area += a;
            } else {
                ocean_area += a;
            }
        }
        let frac = land_area / (land_area + ocean_area + 1e-12);
        println!("[diag] hyps: land={:.1}%", 100.0 * frac);
    }

    // Advect scalar/rheologic fields under rigid plate motion
    // 1) C and th_c_m via existing helper
    crate::continent::advect_c_thc(
        &world.grid,
        &world.v_en,
        dt as f64,
        &mut world.c,
        &mut world.th_c_m,
    );
    // Clamp C to [0,1] defensively
    for v in &mut world.c {
        *v = v.clamp(0.0, 1.0);
    }
    // 2) Age, sediment via generic scalar advection using 3D velocities
    // Build 3D velocities if not already (vel3 above reflects enable_rigid_motion)
    let vel3_for_adv =
        if cfg.enable_rigid_motion { vel3.clone() } else { vec![[0.0; 3]; world.grid.cells] };
    {
        let mut age_next = world.age_myr.clone();
        crate::sl_advect::advect_scalar(
            &world.grid,
            &vel3_for_adv,
            dt as f64,
            &world.age_myr,
            &mut age_next,
        );
        world.age_myr = age_next;
    }
    {
        let mut sed_next = world.sediment_m.clone();
        crate::sl_advect::advect_scalar(
            &world.grid,
            &vel3_for_adv,
            dt as f64,
            &world.sediment_m,
            &mut sed_next,
        );
        world.sediment_m = sed_next;
    }

    // Kinematics/advection (plates velocities already in world.v_en); age increments handled by ridge stage below
    // For MVP, reuse existing world::step_once pieces indirectly through ridge/subduction updates

    // Ridge births and freshening ages
    let mut ages = world.age_myr.clone();
    let _ridge_stats = ridge::apply_ridge(
        &world.grid,
        &world.boundaries,
        &mut ages,
        ridge::RidgeParams { fringe_age_myr: 0.0 },
    );
    // Age growth except ridges set to zero
    for (aw, ar) in world.age_myr.iter_mut().zip(ages.iter()) {
        if *ar == 0.0 {
            *aw = 0.0;
        } else {
            *aw += dt;
        }
    }

    // Build age-only baselines: previous and current
    let n_cells = world.grid.cells;
    let dt_f = dt as f64;
    let mut d0_prev = vec![0.0f32; n_cells];
    let mut d0_now = vec![0.0f32; n_cells];
    for i in 0..n_cells {
        // Approximate previous age: current minus dt unless reset at ridge (ages[i]==0)
        let age_now = world.age_myr[i] as f64;
        let age_was = if ages[i] == 0.0 { 0.0 } else { (age_now - dt_f).max(0.0) };
        let prev = crate::age::depth_from_age_blend(age_was, 2500.0, 350.0, 6300.0, 60.0) as f32;
        let curr = crate::age::depth_from_age_blend(age_now, 2500.0, 350.0, 6300.0, 60.0) as f32;
        d0_prev[i] = prev;
        d0_now[i] = curr;
    }
    // Compose tectonic edits as a single delta added to baseline
    let n = world.grid.cells;
    let base_depth = world.depth_m.clone(); // previous full depth (for dv_struct)
    let th_before = world.th_c_m.clone();
    let mut delta_tect: Vec<f32> = vec![0.0; n];
    let cap_per_step: f32 = 200.0; // global trust region cap (meters)
    if cfg.enable_subduction && do_sub {
        // Adaptive substeps to satisfy CFL ≤ 0.5 for trench normal speeds
        let mut n_sub = cfg.substeps_subduction.max(1);
        // Estimate CFL using max normal speed over convergent edges
        let mut vmax_n: f64 = 0.0;
        for ek in &world.boundaries.edge_kin {
            if ek.class as u8 == 2 {
                vmax_n = vmax_n.max((ek.n_m_per_yr as f64).abs());
            }
        }
        let w_half_m = (cfg.sub_trench_half_width_km as f64).max(1.0) * 1000.0;
        if w_half_m > 0.0 {
            let cfl = (vmax_n * (dt as f64 * 1.0e6)) / w_half_m;
            if cfl > 0.5 {
                let needed = (cfl / 0.5).ceil() as u32;
                let new_sub = needed.clamp(n_sub, 32);
                if new_sub > n_sub {
                    println!(
                        "[substep] subduction: increasing substeps {} → {} (CFL {:.2})",
                        n_sub, new_sub, cfl
                    );
                    n_sub = new_sub;
                }
            }
        }
        let sub_dt = (dt as f64) / (n_sub as f64);
        let cap: f32 = 200.0; // per-substep cap
        let mut cap_count_total = 0u32;
        let mut units_bug_max = 0.0f32;
        for _ in 0..n_sub {
            let mut sub_delta = vec![0.0f32; n];
            let _sub_stats = subduction::compute_subduction_delta(
                &world.grid,
                &world.boundaries,
                &world.plates.plate_id,
                &world.age_myr,
                &world.v_en,
                &mut sub_delta,
                subduction::SubductionParams {
                    tau_conv_m_per_yr: cfg.sub_tau_conv_m_per_yr as f64,
                    trench_half_width_km: 40.0,
                    arc_offset_km: 140.0,
                    arc_half_width_km: 30.0,
                    backarc_width_km: 150.0,
                    trench_deepen_m: 20.0,
                    arc_uplift_m: -12.0,
                    backarc_uplift_m: -10.0,
                    rollback_offset_m: cfg.sub_rollback_offset_m as f64,
                    rollback_rate_km_per_myr: cfg.sub_rollback_rate_km_per_myr as f64,
                    backarc_extension_mode: cfg.sub_backarc_extension_mode,
                    backarc_extension_deepen_m: cfg.sub_backarc_extension_deepen_m,
                    continent_c_min: 0.60,
                },
                Some(&world.c),
            );
            let mut cap_count = 0u32;
            // CFL diagnostic: v*dt / w_half
            let _vmax_t_m_per_yr: f64 = 0.0;
            let mut vmax_n_m_per_yr: f64 = 0.0;
            for ek in &world.boundaries.edge_kin {
                if ek.class as u8 == 2 {
                    vmax_n_m_per_yr = vmax_n_m_per_yr.max((ek.n_m_per_yr as f64).abs());
                }
            }
            // Estimate half-width in meters from cfg
            let w_half_m = (cfg.sub_trench_half_width_km as f64) * 1000.0;
            let cfl =
                if w_half_m > 0.0 { (vmax_n_m_per_yr * (sub_dt * 1.0e6)) / w_half_m } else { 0.0 };
            if cfl > 0.5 && (k % 8 == 0) {
                println!("[cfl] subduction: v_n*dt/w_half = {:.2} > 0.5 (v_n={:.3} m/yr, dt={:.2} Myr, w_half={:.0} m)", cfl, vmax_n_m_per_yr, sub_dt, w_half_m);
            }
            for i in 0..n {
                let raw = sub_delta[i] * (sub_dt as f32);
                if raw.abs() > 1000.0 {
                    units_bug_max = units_bug_max.max(raw.abs());
                }
                let clamped = raw.clamp(-cap, cap);
                if raw.abs() > cap {
                    cap_count += 1;
                }
                delta_tect[i] += clamped;
            }
            cap_count_total += cap_count;
        }
        if cap_count_total > 0 {
            println!(
                "[cap] subduction: {} hits across substeps capped to ±{:.0} m/substep",
                cap_count_total, 200.0
            );
        }
        if units_bug_max > 0.0 {
            println!("[UNITS_BUG] subduction Δz max {:.1} m (>1000 m)", units_bug_max);
        }
    }
    // Continental rifting (delta-only) before transforms
    {
        let p_rift = crate::rifting::RiftingParams {
            c_rift_min: 0.60,
            v_open_min_m_per_yr: 0.015, // 15 mm/yr
            w_core_km: 60.0,
            w_taper_km: 200.0,
            k_thin: 0.0025,
            alpha_subs: 0.33,
            ocean_thresh: 0.15,
            k_c_oceanize: 0.010,
            reset_age_on_core: true,
            enable_shoulder: true,
            w_bulge_km: 80.0,
            beta_shoulder: 0.10,
            couple_flexure: false,
            thc_min_m: 25_000.0,
            thc_max_m: 65_000.0,
        };
        let mut rift_tmp = vec![0.0f32; n];
        let _ = crate::rifting::apply_rifting(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &mut world.c,
            &mut world.th_c_m,
            &mut world.age_myr,
            &mut rift_tmp,
            &world.area_m2,
            &p_rift,
            dt as f64,
        );
        let mut units_bug_max = 0.0f32;
        for i in 0..n {
            let raw = rift_tmp[i];
            if raw.abs() > 1000.0 {
                units_bug_max = units_bug_max.max(raw.abs());
            }
            let add = raw.clamp(-cap_per_step, cap_per_step);
            delta_tect[i] += add;
        }
        if units_bug_max > 0.0 {
            println!("[UNITS_BUG] rifting Δz max {:.1} m (>1000 m)", units_bug_max);
        }
    }
    // Transforms after rifting
    if do_trf {
        // Adaptive substeps to satisfy CFL ≤ 0.5 for tangential speeds
        let mut n_sub = cfg.substeps_transforms.max(1);
        let mut vmax_t: f64 = 0.0;
        for ek in &world.boundaries.edge_kin {
            if ek.class as u8 == 3 {
                vmax_t = vmax_t.max((ek.t_m_per_yr as f64).abs());
            }
        }
        let w_half_m = 1000.0 * 10.0; // match basin_half_width_km below
        let cfl = if w_half_m > 0.0 { (vmax_t * (dt as f64 * 1.0e6)) / w_half_m } else { 0.0 };
        if cfl > 0.5 {
            let needed = (cfl / 0.5).ceil() as u32;
            let new_sub = needed.clamp(n_sub, 32);
            if new_sub > n_sub {
                println!(
                    "[substep] transforms: increasing substeps {} → {} (CFL {:.2})",
                    n_sub, new_sub, cfl
                );
                n_sub = new_sub;
            }
        }
        let sub_dt = (dt as f64) / (n_sub as f64);
        let cap: f32 = 200.0;
        let mut cap_hits = 0u32;
        let mut units_bug_max = 0.0f32;
        for _ in 0..n_sub {
            let mut tr_tmp = vec![0.0f32; n];
            let _ = transforms::apply_transforms(
                &world.grid,
                &world.boundaries,
                &world.plates.plate_id,
                &world.v_en,
                &mut tr_tmp,
                transforms::TransformParams {
                    tau_open_m_per_yr: 0.005,
                    min_tangential_m_per_yr: 0.015,
                    max_normal_m_per_yr: 0.001,
                    basin_half_width_km: 10.0,
                    ridge_like_uplift_m: -3.0,
                    basin_deepen_m: 6.0,
                },
                sub_dt,
            );
            // CFL diagnostic for transforms using tangential speed over half-width
            let mut vmax_t_m_per_yr: f64 = 0.0;
            for ek in &world.boundaries.edge_kin {
                if ek.class as u8 == 3 {
                    vmax_t_m_per_yr = vmax_t_m_per_yr.max((ek.t_m_per_yr as f64).abs());
                }
            }
            let w_half_m = 1000.0 * 10.0; // basin_half_width_km default above
            let cfl = (vmax_t_m_per_yr * (sub_dt * 1.0e6)) / w_half_m;
            if cfl > 0.5 && (k % 8 == 0) {
                println!("[cfl] transforms: v_t*dt/w_half = {:.2} > 0.5 (v_t={:.3} m/yr, dt={:.2} Myr, w_half={:.0} m)", cfl, vmax_t_m_per_yr, sub_dt, w_half_m);
            }
            for i in 0..n {
                let raw = tr_tmp[i];
                if raw.abs() > 1000.0 {
                    units_bug_max = units_bug_max.max(raw.abs());
                }
                let clamped = raw.clamp(-cap, cap);
                if raw.abs() > cap {
                    cap_hits += 1;
                }
                delta_tect[i] += clamped;
            }
        }
        if cap_hits > 0 {
            println!(
                "[cap] transforms: {} hits across substeps capped to ±{:.0} m/substep",
                cap_hits, 200.0
            );
        }
        if units_bug_max > 0.0 {
            println!("[UNITS_BUG] transforms Δz max {:.1} m (>1000 m)", units_bug_max);
        }
    }
    // Continents buoyancy: drive persistent delta_buoy_m toward target amplitude, compose as delta
    {
        // Densities
        let rho_m = 3300.0f32;
        let rho_c = 2850.0f32;
        let rho_w = 1000.0f32;
        let rho_a = 1.2f32;
        let mut dstep_max = 0.0f32;
        let mut amp_max = 0.0f32;
        let step_cap = cap_per_step; // reuse global per-step cap
        if world.delta_buoy_m.len() != n { world.delta_buoy_m.resize(n, 0.0); }
        for (i, dt) in delta_tect.iter_mut().enumerate().take(n) {
            let c = world.c[i].clamp(0.0, 1.0);
            let th = world.th_c_m[i].clamp(25_000.0, 65_000.0);
            let elev_prev = world.sea.eta_m - world.depth_m[i];
            let rho_top = if elev_prev > 0.0 { rho_a } else { rho_w };
            let frac = (rho_m - rho_c) / (rho_m - rho_top);
            let target = -c * frac * (th - 35_000.0);
            let old = world.delta_buoy_m[i];
            let raw = target - old;
            let step = raw.clamp(-step_cap, step_cap);
            let new_amp = old + step;
            world.delta_buoy_m[i] = new_amp;
            dstep_max = dstep_max.max(step.abs());
            amp_max = amp_max.max(new_amp.abs());
            *dt += step; // apply only the change this step
        }
        println!("[buoyancy] dZ_step max={:.1} m | amp max={:.0} m", dstep_max, amp_max);
        world.epoch_continents = world.epoch_continents.wrapping_add(1);
    }
    // Accretion and Orogeny deltas after continents uplift
    {
        // Reuse last subduction masks by recomputing quickly if needed
        let mut sub_tmp = vec![0.0f32; n];
        let sub = subduction::apply_subduction(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            &mut sub_tmp,
            subduction::SubductionParams {
                tau_conv_m_per_yr: cfg.sub_tau_conv_m_per_yr as f64,
                trench_half_width_km: cfg.sub_trench_half_width_km as f64,
                arc_offset_km: cfg.sub_arc_offset_km as f64,
                arc_half_width_km: cfg.sub_arc_half_width_km as f64,
                backarc_width_km: cfg.sub_backarc_width_km as f64,
                trench_deepen_m: cfg.sub_trench_deepen_m,
                arc_uplift_m: cfg.sub_arc_uplift_m,
                backarc_uplift_m: cfg.sub_backarc_uplift_m,
                rollback_offset_m: cfg.sub_rollback_offset_m as f64,
                rollback_rate_km_per_myr: cfg.sub_rollback_rate_km_per_myr as f64,
                backarc_extension_mode: cfg.sub_backarc_extension_mode,
                backarc_extension_deepen_m: cfg.sub_backarc_extension_deepen_m,
                continent_c_min: cfg.sub_continent_c_min,
            },
            Some(&world.c),
        );
        // O–C accretion
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
            &sub.masks,
            &world.boundaries,
            &vel3,
            &mut world.c,
            &mut world.th_c_m,
            &mut acc_tmp,
            &world.area_m2,
            &p_acc,
            dt as f64,
        );
        let mut units_bug_max_acc = 0.0f32;
        for i in 0..n {
            let raw = acc_tmp[i];
            if raw.abs() > 1000.0 {
                units_bug_max_acc = units_bug_max_acc.max(raw.abs());
            }
            delta_tect[i] += raw.clamp(-cap_per_step, cap_per_step);
        }
        if units_bug_max_acc > 0.0 {
            println!("[UNITS_BUG] accretion Δz max {:.1} m (>1000 m)", units_bug_max_acc);
        }
        // C–C orogeny
        let p_orog = crate::orogeny::OrogenyParams {
            c_min: 0.60,
            w_core_km: 100.0,
            w_taper_km: 200.0,
            k_thick: 0.002,
            beta_uplift: 0.30,
            gamma_obliquity: 1.5,
            couple_flexure: false,
        };
        let mut orog_tmp = vec![0.0f32; n];
        let _ = crate::orogeny::apply_cc_orogeny(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.c,
            &world.area_m2,
            &mut world.th_c_m,
            &mut orog_tmp,
            &p_orog,
            dt as f64,
        );
        let mut units_bug_max_orog = 0.0f32;
        for i in 0..n {
            let raw = orog_tmp[i];
            if raw.abs() > 1000.0 {
                units_bug_max_orog = units_bug_max_orog.max(raw.abs());
            }
            delta_tect[i] += raw.clamp(-300.0, 300.0);
        }
        if units_bug_max_orog > 0.0 {
            println!("[UNITS_BUG] orogeny Δz max {:.1} m (>1000 m)", units_bug_max_orog);
        }

        // Composite per-cell cap on the sum across operators
        let mut composite_caps = 0u32;
        for d in &mut delta_tect {
            if d.abs() > cap_per_step {
                *d = d.clamp(-cap_per_step, cap_per_step);
                composite_caps += 1;
            }
        }
        if composite_caps > 0 {
            println!(
                "[cap] composite per-cell: clamped {} cells (cap={:.1} m)",
                composite_caps, cap_per_step
            );
        }

        // Global trust region: scale deltas if any exceed per-step hard limit
        let mut max_abs = 0.0f32;
        for &d in &delta_tect {
            max_abs = max_abs.max(d.abs());
        }
        if max_abs > cap_per_step {
            let scale = cap_per_step / max_abs;
            for d in &mut delta_tect {
                *d *= scale;
            }
            println!("[trust] scaled tectonic deltas by {:.3} (max {:.1} m)", scale, max_abs);
        }
        // Enforce per-step cap on crustal thickness changes
        let mut th_caps = 0u32;
        let mut th_units_bug = 0.0f32;
        for (i, _) in th_before.iter().enumerate().take(n) {
            let dth = world.th_c_m[i] - th_before[i];
            if dth.abs() > 1000.0 {
                th_units_bug = th_units_bug.max(dth.abs());
            }
            if dth > cap_per_step {
                world.th_c_m[i] = th_before[i] + cap_per_step;
                th_caps += 1;
            }
            if dth < -cap_per_step {
                world.th_c_m[i] = th_before[i] - cap_per_step;
                th_caps += 1;
            }
        }
        if th_caps > 0 {
            println!("[cap] th_c per-step: clamped {} cells (cap={:.1} m)", th_caps, cap_per_step);
        }
        if th_units_bug > 0.0 {
            println!("[UNITS_BUG] Δth_c max {:.1} m (>1000 m)", th_units_bug);
        }
    }
    // Preserve previous structural component relative to age baseline, then add incremental deltas
    for (i, _) in th_before.iter().enumerate().take(n) {
        let struct_prev = base_depth[i] - d0_prev[i];
        world.depth_m[i] = d0_now[i] + struct_prev + delta_tect[i];
    }

    // Flexure: experimental GPU V-cycle or Winkler fallback
    if cfg.enable_flexure && do_flx {
        // Update Te field from age/C/th_c before assembling load
        crate::flexure::compute_te_field(&world.age_myr, &world.c, &world.th_c_m, &mut world.te_m);
        let lp = flexure_loads::LoadParams {
            rho_w: 1030.0,
            rho_c: 2850.0,
            g: 9.81,
            sea_level_m: *surf.eta_m,
        };
        // Use compaction-aware sediment load if available
        let mut f_load = flexure_loads::assemble_load_with_sediments(
            &world.grid,
            &world.depth_m,
            &world.sediment_m,
            &lp,
        );
        if cfg.subtract_mean_load {
            let mut sum: f64 = 0.0;
            for &q in &f_load {
                sum += q as f64;
            }
            let mean = (sum / (f_load.len().max(1) as f64)) as f32;
            for q in &mut f_load {
                *q -= mean;
            }
        }
        let n = world.grid.cells;
        if cfg.use_gpu_flexure {
            // Best-effort GPU path: use engine::gpu context and flexure_gpu V-cycle
            // Note: This creates a transient context per call; future work should reuse.
            let fut = async {
                let ctx = crate::gpu::persistent();
                let mut flex = crate::flexure_gpu::FlexGpu::new(ctx);
                // Tile approximation over a rectangular patch (no halos in interior), dx ≈ mean cell scale
                let mut area_mean: f64 = 0.0;
                for a in &world.area_m2 {
                    area_mean += *a as f64;
                }
                let cell_scale_m = (area_mean.max(1.0) / (n as f64)).sqrt() as f32;
                // Choose a near-square layout to reduce aspect ratio artifacts
                let w_guess = (n as f64).sqrt().floor() as u32;
                let width = w_guess.clamp(64, 4096).max(1);
                let height = ((n as u32) + width - 1) / width;
                let dims =
                    crate::flexure_gpu::TileDims { width, height, halo: 2, dx: cell_scale_m };
                let w_tot = dims.width_total() as usize;
                let h_tot = dims.height_total() as usize;
                let n_tot = w_tot * h_tot;
                let mut w_tex = crate::flexure_gpu::GpuTex::new_storage(ctx, n_tot, "flex.w");
                let f_tex = crate::flexure_gpu::GpuTex::new_storage(ctx, n_tot, "flex.f");
                // Stage f into interior region row-by-row with halos
                let mut f_stage = vec![0.0f32; n_tot];
                let halo = dims.halo as usize;
                let interior_w = dims.width as usize;
                let interior_h = dims.height as usize;
                let mut src_off = 0usize;
                for iy in 0..interior_h {
                    let dst_row = iy + halo;
                    let dst_start = dst_row * w_tot + halo;
                    let n_copy = interior_w.min(n - src_off);
                    if n_copy > 0 {
                        f_stage[dst_start..dst_start + n_copy]
                            .copy_from_slice(&f_load[src_off..src_off + n_copy]);
                        src_off += n_copy;
                    }
                    if src_off >= n {
                        break;
                    }
                }
                ctx.queue.write_buffer(&f_tex.buf, 0, bytemuck::cast_slice(&f_stage));
                let mut scratch = crate::flexure_gpu::FlexScratch::new(ctx, n_tot);
                // Use area-weighted mean Te to compute a single D for this MVP
                let mut sum_a = 0.0f64;
                let mut sum_te = 0.0f64;
                for i in 0..n {
                    let a = world.area_m2[i] as f64;
                    sum_a += a;
                    sum_te += a * world.te_m[i] as f64;
                }
                let te_mean = if sum_a > 0.0 { (sum_te / sum_a) as f32 } else { 25_000.0 };
                let e_pa = 70.0e9f32;
                let nu = 0.25f32;
                let d = e_pa * te_mean.powi(3) / (12.0 * (1.0 - nu * nu));
                let p = crate::flexure_gpu::FlexParams {
                    d,
                    k: 0.0,
                    wj_omega: cfg.gpu_wj_omega.clamp(0.4, 0.95),
                    nu1: 1,
                    nu2: 1,
                    levels: cfg.gpu_flex_levels.max(1),
                };
                // Run requested number of V-cycles (reduced to 1 for stability)
                for _ in 0..1u32 {
                    let _stats = flex.v_cycle(ctx, dims, &mut w_tex, &f_tex, &mut scratch, &p);
                }
                // Read back w and extract interior row
                let w_all =
                    crate::flexure_gpu::blocking_read(&ctx.device, &ctx.queue, &w_tex.buf, n_tot);
                let mut w_out = vec![0.0f32; n];
                let mut dst_off = 0usize;
                for iy in 0..interior_h {
                    let src_row = iy + halo;
                    let src_start = src_row * w_tot + halo;
                    let n_copy = interior_w.min(n - dst_off);
                    if n_copy > 0 {
                        w_out[dst_off..dst_off + n_copy]
                            .copy_from_slice(&w_all[src_start..src_start + n_copy]);
                        dst_off += n_copy;
                    }
                    if dst_off >= n {
                        break;
                    }
                }
                w_out
            };
            let w_host: Vec<f32> = pollster::block_on(fut);
            for (d, wi) in world.depth_m.iter_mut().zip(w_host.iter()) {
                *d = (*d + *wi).clamp(-8000.0, 8000.0);
            }
        } else {
            let k = 3.0e8f32; // Winkler fallback
            for (d, &q) in world.depth_m.iter_mut().zip(f_load.iter()) {
                let w = q / k;
                *d = (*d + w).clamp(-8000.0, 8000.0);
            }
        }
    }
    // NOTE(Science): Winkler foundation ignores lithospheric flexural coupling length scales.
    // It is useful as a placeholder but cannot reproduce forebulge/trough patterns.

    // Erosion/diffusion
    let do_surf = k % (cfg.cadence_surf_every.max(1) as u64) == 0;
    if cfg.enable_erosion && do_surf {
        let sp = crate::surface::SurfaceParams {
            k_stream: cfg.surf_k_stream,
            m_exp: cfg.surf_m_exp,
            n_exp: cfg.surf_n_exp,
            k_diff: cfg.surf_k_diff,
            k_tr: cfg.surf_k_tr,
            p_exp: cfg.surf_p_exp,
            q_exp: cfg.surf_q_exp,
            rho_sed: cfg.surf_rho_sed,
            min_slope: cfg.surf_min_slope,
            subcycles: cfg.surf_subcycles.max(1),
            couple_flexure: cfg.surf_couple_flexure,
        };
        let _stats = crate::surface::apply_surface_processes(
            &world.grid,
            &world.c,
            &mut world.depth_m,
            &mut world.sediment_m,
            &world.area_m2,
            &sp,
            cfg.dt_myr as f64,
        );
    }

    // Solve eta to hit target land fraction unless frozen; elevation is independent of eta
    // Eta solve moved after writing elevation to reduce coastline speckle

    if !cfg.freeze_eta && do_sea {
        // Maintain constant ocean volume policy (isostasy) solving directly for η near current value
        if world.sea_level_ref.is_none() {
            // Initialize eta by target land fraction via hypsometry (avoids 0%/100% starts)
            let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
            let target_ocean = (1.0 - cfg.target_land_frac.clamp(0.0, 1.0)) as f64;
            let (eta0, ocean_area_m2, volume_m3) =
                crate::isostasy::solve_eta_for_ocean_area_fraction_hypso(
                    &world.depth_m,
                    &world.area_m2,
                    target_ocean,
                );
            *surf.eta_m = eta0 as f32;
            world.sea_level_ref = Some(crate::world::SeaLevelRef { volume_m3, ocean_area_m2 });
            let frac = if total_area > 0.0 { ocean_area_m2 / total_area } else { 0.0 };
            println!("[isostasy] init: ocean_frac={:.3} eta0={:.0} m", frac, eta0);
        }
        if let Some(r) = world.sea_level_ref {
            // Sanitize inputs
            for d in &mut world.depth_m {
                if !d.is_finite() {
                    *d = 0.0;
                }
                *d = d.clamp(-9000.0, 11000.0);
            }
            // Hypsometry-based exact eta solve with guardrails
            let eta_new = crate::isostasy::solve_eta_hypsometry(
                &world.depth_m,
                &world.area_m2,
                r.volume_m3,
                *surf.eta_m as f64,
            );
            // Δη clamp by structural volume change estimate
            let mut dv_pos = 0.0f64; // accommodation created
            let mut dv_neg = 0.0f64; // infilling
            for (a_d0, d1) in world.area_m2.iter().zip(base_depth.iter()).zip(world.depth_m.iter())
            {
                let (a, d0) = a_d0;
                let dd = (*d1 as f64) - (*d0 as f64);
                if dd > 0.0 {
                    dv_pos += (*a as f64) * dd;
                } else {
                    dv_neg += (*a as f64) * (-dd);
                }
            }
            let dv_struct = dv_pos - dv_neg;
            // Approx ocean area ~ half total if unknown; safer: compute at prev eta
            let mut a_ocean = 0.0f64;
            for (a, d) in world.area_m2.iter().zip(world.depth_m.iter()) {
                if (*d as f64 - *surf.eta_m as f64) > 0.0 {
                    a_ocean += *a as f64;
                }
            }
            let a_ref = a_ocean.max(1.0e6);
            let d_eta_allow = (dv_struct / a_ref).clamp(-200.0, 200.0);
            let eta_clamped =
                eta_new.clamp((*surf.eta_m as f64) - 200.0, (*surf.eta_m as f64) + 200.0);
            let eta_final = if (eta_clamped - *surf.eta_m as f64).abs() <= d_eta_allow.abs() {
                eta_clamped
            } else {
                *surf.eta_m as f64 + d_eta_allow
            };
            *surf.eta_m = eta_final as f32;
        }
    }

    // Write surface elevation after eta solve: z = eta - depth
    for (e, &d) in surf.elev_m.iter_mut().zip(world.depth_m.iter()) {
        let z = (*surf.eta_m as f64 - d as f64) as f32;
        *e = if z.is_finite() { z } else { 0.0 };
    }

    // Sanitize
    crate::util::sanitize_elevation(surf.elev_m);

    // Optional: mass/volume budgets and per-step deltas
    if cfg.log_mass_budget {
        let rho_cont = 2900.0f64;
        // Effective sediment density with compaction for budgets
        let mut m_cont: f64 = 0.0;
        let mut m_sed: f64 = 0.0;
        let mut vol_ocean: f64 = 0.0;
        let mut land_area: f64 = 0.0;
        let mut tot_area: f64 = 0.0;
        for i in 0..world.grid.cells {
            let a = world.area_m2[i] as f64;
            let d = world.depth_m[i] as f64;
            tot_area += a;
            m_cont += a * (world.th_c_m.get(i).copied().unwrap_or(0.0) as f64) * rho_cont;
            let sed_thk = world.sediment_m.get(i).copied().unwrap_or(0.0);
            let rho_eff = crate::flexure_loads::sediment_rho_eff(sed_thk) as f64;
            m_sed += a * (sed_thk as f64) * rho_eff;
            // Ocean thickness = max(depth - eta, 0)
            let h_ocean = (d - *surf.eta_m as f64).max(0.0);
            vol_ocean += a * h_ocean;
            if (*surf.eta_m as f64 - d) > 0.0 {
                land_area += a;
            }
        }
        let land_frac = if tot_area > 0.0 { land_area / tot_area } else { 0.0 };
        let d_cont = m_cont - world.last_mass_cont_kg;
        let d_sed = m_sed - world.last_mass_sed_kg;
        let d_vol = vol_ocean - world.last_ocean_vol_m3;
        // Reference ocean volume if available
        let vol_err = if let Some(r) = world.sea_level_ref { vol_ocean - r.volume_m3 } else { 0.0 };
        println!(
            "[budget] land={:.1}% M_cont={:.3e} (d={:+.1e}) M_sed={:.3e} (d={:+.1e}) V_ocean={:.3e} (d={:+.1e}) err_ref={:+.1e}",
            100.0 * land_frac,
            m_cont,
            d_cont,
            m_sed,
            d_sed,
            vol_ocean,
            d_vol,
            vol_err
        );
        world.last_mass_cont_kg = m_cont;
        world.last_mass_sed_kg = m_sed;
        world.last_ocean_vol_m3 = vol_ocean;
    }
    // Advance clock for viewer pipeline path
    world.clock.t_myr += dt as f64;
    world.clock.step_idx = world.clock.step_idx.saturating_add(1);
}
