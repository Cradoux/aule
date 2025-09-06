//! Centralized physics pipeline (T-632): run full step sequence and expose a single surface elevation.
//!
//! This viewer-facing pipeline contrasts with `world::step_once` in one key way:
//! - Here we compute and output `eta` (sea level) but do not bake the offset back into `world.depth_m`.
//!   Callers must render elevation as `z = -depth - eta`. This reduces coastline flicker and keeps
//!   `depth_m` as a separable tectonic + age baseline.
//! - Subduction/transform/rifting/orogeny still edit `depth_m` additively on top of age-derived bathy.
//!   If these are applied at high cadence without normalization, long-run drift of `depth_m` can occur.

use crate::{boundaries, config::PipelineCfg, flexure_loads, ridge, subduction, transforms, world::World};
use std::thread;


/// Borrowed views of surface fields required by the pipeline.
pub struct SurfaceFields<'a> {
    /// Solved surface elevation (meters), vertex-major, same order as GPU uploads.
    pub elevation_m: &'a mut [f32],
    /// Sea level eta (meters) relative to geoid; output of sea-level solve.
    pub eta_m: &'a mut f32,
}

/// Invariants & guards (TST-2):
/// - No negative crustal thickness (`th_c_m >= 0`)
/// - Slopes within reasonable bounds when provided elevation (unitless slope ≤ 10)
///
/// In debug builds, violations panic via debug_assert. In release, they are logged with counts.
pub(crate) fn check_invariants(world: &World, elev_opt: Option<&[f32]>) {
    // 1) No negative crustal thickness
    let mut neg_th = 0u32;
    for &t in &world.th_c_m {
        if t < 0.0 {
            neg_th += 1;
        }
    }
    if neg_th > 0 {
        println!("[invar][th_c] negative_count={}", neg_th);
    }
    debug_assert_eq!(neg_th, 0, "Negative crustal thickness encountered");

    // 2) Slope bound check if elevation provided
    if let Some(elev) = elev_opt {
        let n = world.grid.cells.min(elev.len());
        let r = 6_371_000.0_f64;
        let mut over_slope = 0u32;
        for i in 0..n {
            let lat_i = world.grid.latlon[i][0] as f64;
            let lon_i = world.grid.latlon[i][1] as f64;
            let zi = elev[i] as f64;
            if !zi.is_finite() {
                continue;
            }
            if zi.abs() > 20_000.0 {
                continue;
            }
            let cos_lat = lat_i.cos();
            let mut s_max = 0.0_f64;
            for &nj in &world.grid.n1[i] {
                let j = nj as usize;
                if j >= n {
                    continue;
                }
                let lat_j = world.grid.latlon[j][0] as f64;
                let lon_j = world.grid.latlon[j][1] as f64;
                let dx = (lon_j - lon_i) * cos_lat * r; // east (m)
                let dy = (lat_j - lat_i) * r; // north (m)
                let dist = (dx * dx + dy * dy).sqrt();
                if dist > 0.0 {
                    let zj = elev[j] as f64;
                    if !zj.is_finite() {
                        continue;
                    }
                    if zj.abs() > 20_000.0 {
                        continue;
                    }
                    let s = (zj - zi).abs() / dist;
                    if s > s_max {
                        s_max = s;
                    }
                }
            }
            if s_max > 10.0 {
                over_slope += 1;
            }
        }
        if over_slope > 0 {
            println!("[invar][slope] cells_over_10={}", over_slope);
        }
        // Keep non-fatal: log only; slope spikes may occur with outlier elevations
    }
}

/// DEPRECATED: Use unified_pipeline::UnifiedPipeline instead.
/// This function is kept temporarily for tests but will be removed.
#[deprecated(note = "Use unified_pipeline::UnifiedPipeline::step instead")]
pub fn step_full(world: &mut World, surf: SurfaceFields, cfg: PipelineCfg) {
    // Performance timers (ms)
    let mut ms_boundaries = 0.0f64;
    let mut ms_kinematics = 0.0f64;
    let mut ms_ridge = 0.0f64;
    let mut ms_subduction = 0.0f64;
    let mut ms_transforms = 0.0f64;
    let mut ms_buoyancy = 0.0f64;
    let mut ms_accretion_orogeny = 0.0f64;
    let mut ms_flexure = 0.0f64;
    let mut ms_surface = 0.0f64;
    let mut ms_eta = 0.0f64;
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
        let t_kin0 = std::time::Instant::now();
        // One-shot remap: snapshot fields, substep plate motion only, then remap once
        let dt_years = (cfg.dt_myr as f64) * 1.0e6;
        let c_src = world.c.clone();
        let th_src = world.th_c_m.clone();
        let age_src = world.age_myr.clone();
        
        // Debug: Check if continents exist before advection (for mass conservation)
        let c_sum_before = c_src.iter().sum::<f32>();

        let theta_max = world.plates.max_abs_omega_rad_yr() * dt_years;
        let cell_angle = world.grid.mean_cell_angle_rad();
        let theta_limit = 1.0 * cell_angle; // allow up to one cell crossing per substep
        let nsub = ((theta_max / theta_limit).ceil() as u32).clamp(1, 8);
        for _ in 0..nsub {
            world.plates.advance_rigid(&world.grid, dt_years / (nsub as f64));
        }

        // One-shot remap from snapshots
        crate::continent::rigid_advect_c_thc_from(
            &world.grid,
            &world.plates,
            dt_years,
            &c_src,
            &th_src,
            &mut world.c,
            &mut world.th_c_m,
        );
        
        // Check for mass conservation and apply correction if needed
        let c_sum_after = world.c.iter().sum::<f32>();
        
        // Debug: Check if plate motion is working (always show to track force balance)
        let max_omega = world.plates.omega_rad_yr.iter().fold(0.0f32, |a, &b| a.max(b.abs()));
        let max_velocity = world.v_en.iter().fold(0.0f32, |a, v| a.max((v[0]*v[0] + v[1]*v[1]).sqrt()));
        println!("[rigid_motion] BEFORE_FB: max_omega={:.6} rad/yr, max_velocity={:.6} m/yr", max_omega, max_velocity);
        
        // CRITICAL FIX: Enforce mass conservation to prevent continental loss
        let mass_loss = c_sum_before - c_sum_after;
        if mass_loss.abs() > 0.1 {
            // Silently correct mass loss - this is normal due to interpolation
            // Redistribute lost mass proportionally to existing continents
            let total_existing = c_sum_after;
            if total_existing > 0.0 {
                let correction_factor = c_sum_before / total_existing;
                for c_val in &mut world.c {
                    if *c_val > 0.0 {
                        *c_val *= correction_factor;
                        *c_val = c_val.clamp(0.0, 1.0); // Keep in valid range
                    }
                }
                // Mass conservation correction applied
            }
        }
        crate::age::rigid_advect_age_from(
            &world.grid,
            &world.plates,
            dt_years,
            &age_src,
            &mut world.age_myr,
        );
        crate::age::increment_oceanic_age(&world.c, &mut world.age_myr, cfg.dt_myr as f64);
        // Keep ridges wet: enforce ocean only where very young seafloor AND previously oceanic
        // CRITICAL FIX: Protect continents from ridge destruction
        for i in 0..world.grid.cells {
            if world.age_myr[i] < 2.0 && c_src.get(i).copied().unwrap_or(0.0) < 0.1 && world.c[i] < 0.1 {
                // Only destroy oceanic areas (both before AND after advection must be oceanic)
                world.c[i] = 0.0;
                world.th_c_m[i] = 0.0;
            }
        }
        
        // Note: Continental uplift will be applied later after all geological processes

        // Rebuild velocities and kinds after motion
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
        vel3 = crate::plates::velocity_field_m_per_yr(
            &world.grid,
            &world.plates,
            &world.plates.plate_id,
        );
        world.plates.kind =
            crate::plates::derive_kinds(&world.grid, &world.plates.plate_id, &world.c, 0.35, 0.5);

        // Conservation logs
        let mut area_before = 0.0_f64;
        let mut area_after = 0.0_f64;
        let mut thc_before = 0.0_f64;
        let mut thc_after = 0.0_f64;
        for i in 0..world.grid.cells.min(world.area_m2.len()) {
            let a = world.area_m2[i] as f64;
            let cb = c_src.get(i).copied().unwrap_or(0.0) as f64;
            let ca = world.c.get(i).copied().unwrap_or(0.0) as f64;
            let tb = th_src.get(i).copied().unwrap_or(0.0) as f64;
            let ta = world.th_c_m.get(i).copied().unwrap_or(0.0) as f64;
            area_before += a * cb;
            area_after += a * ca;
            thc_before += a * tb * cb;
            thc_after += a * ta * ca;
        }
        let pct = |d: f64, b: f64| if b.abs() > 0.0 { 100.0 * d / b } else { 0.0 };
        println!(
            "[rigid][cons] C·A d={:+.3}% th_c·A d={:+.3}%",
            pct(area_after - area_before, area_before),
            pct(thc_after - thc_before, thc_before)
        );

        let max_theta_deg =
            world.plates.max_abs_omega_rad_yr() * dt_years * 180.0 / std::f64::consts::PI;
        println!(
            "[rigid] substeps={} max Δθ={:.3}° mean_cell≈{:.3}°",
            nsub,
            max_theta_deg,
            cell_angle * 180.0 / std::f64::consts::PI
        );
        ms_kinematics += t_kin0.elapsed().as_secs_f64() * 1000.0;
    } else {
        world.v_en.clone_from(&world.plates.vel_en);
        vel3.fill([0.0, 0.0, 0.0]);
    }

    // Boundaries classification (after rigid motion and advection if enabled)
    let _span_bound = tracing::info_span!("boundaries").entered();
    let t_b0 = std::time::Instant::now();
    world.boundaries =
        boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
    ms_boundaries += t_b0.elapsed().as_secs_f64() * 1000.0;
    
    // Debug: Log boundary statistics to check if they're changing
    let stats = &world.boundaries.stats;
    println!("[boundaries] ridge={} subd={} trans={} | total_edges={}", 
             stats.divergent, stats.convergent, stats.transform, world.boundaries.edges.len());

    // CRITICAL FIX: Apply force balance to update plate rotation rates
    // This was missing and is why landmasses don't follow plates!
    let t_fb0 = std::time::Instant::now();
    let fb_params = crate::force_balance::FbParams {
        gain: if cfg.fb_gain == 0.0 || cfg.fb_gain < 1.0e-9 { 1.0e-8 } else { cfg.fb_gain }, // Increase gain to restore motion
        damp_per_myr: if cfg.fb_damp_per_myr == 0.0 { 0.2 } else { cfg.fb_damp_per_myr },
        k_conv: if cfg.fb_k_conv == 0.0 { 1.0 } else { cfg.fb_k_conv },
        k_div: if cfg.fb_k_div == 0.0 { 0.5 } else { cfg.fb_k_div },
        k_trans: if cfg.fb_k_trans == 0.0 { 0.1 } else { cfg.fb_k_trans },
        max_domega: 1.0e-7, // Allow reasonable omega changes per step  
        max_omega: 1.0e-6, // Allow reasonable maximum omega values for visible motion
    };
    // Debug: Show force balance parameters to diagnose the issue
    println!("[force_balance] PARAMS: gain={:.2e}, damp={:.3}, k_conv={:.3}, k_div={:.3}, k_trans={:.3}, max_domega={:.2e}, max_omega={:.2e}", 
             fb_params.gain, fb_params.damp_per_myr, fb_params.k_conv, fb_params.k_div, fb_params.k_trans, fb_params.max_domega, fb_params.max_omega);
    // Store plate_id to avoid borrow checker issues
    let plate_id = world.plates.plate_id.clone();
    crate::force_balance::apply_force_balance(
        &world.grid,
        &world.boundaries,
        &plate_id,
        &mut world.plates,
        &world.grid.area,
        cfg.dt_myr as f64,
        fb_params,
    );
    // Update velocity field after force balance changes omega_rad_yr
    world.v_en = crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
    let ms_force_balance = t_fb0.elapsed().as_secs_f64() * 1000.0;
    
    // Debug: Check if force balance updated omega values
    let max_omega_after = world.plates.omega_rad_yr.iter().fold(0.0f32, |a, &b| a.max(b.abs()));
    let max_velocity_after = world.v_en.iter().fold(0.0f32, |a, v| a.max((v[0]*v[0] + v[1]*v[1]).sqrt()));
    println!("[force_balance] AFTER_FB: max_omega={:.6} rad/yr, max_velocity={:.6} m/yr", max_omega_after, max_velocity_after);
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
    // DL-1: Every 1 Myr, emit boundary composition and cap metrics
    if true {
        // Boundary length by class using edge counts as proxy
        let mut n_r = 0u32;
        let mut n_s = 0u32;
        let mut n_t = 0u32;
        for &(_u, _v, cls) in &world.boundaries.edges {
            match cls {
                1 => n_r += 1,
                2 => n_s += 1,
                3 => n_t += 1,
                _ => {}
            }
        }
        let tot = (n_r + n_s + n_t).max(1) as f64;
        let pr = 100.0 * (n_r as f64) / tot;
        let ps = 100.0 * (n_s as f64) / tot;
        let pt = 100.0 * (n_t as f64) / tot;
        // Oceanized cores count (approx: C crossing ocean_thresh)
        let mut oceanized = 0u32;
        for &c in &world.c {
            if c < 0.15 {
                oceanized += 1;
            }
        }
        // Cap metrics: rely on RL-2 softening counters we emit within rifting/orogeny sections
        println!(
            "[diag] boundaries: ridge={:.0}% subd={:.0}% trans={:.0}% | oceanized cores={}",
            pr, ps, pt, oceanized
        );
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
    // Skip legacy semi-Lagrangian advection for C/th_c_m and age when rigid motion is enabled
    if !cfg.enable_rigid_motion {
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
        // 2) Age via generic scalar advection (non-rigid mode only)
        let vel3_for_adv = vec![[0.0; 3]; world.grid.cells];
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
    // Sediment advection still uses the current 3D field
    let vel3_for_adv =
        if cfg.enable_rigid_motion { vel3.clone() } else { vec![[0.0; 3]; world.grid.cells] };
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
    let t_r0 = std::time::Instant::now();
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
    ms_ridge += t_r0.elapsed().as_secs_f64() * 1000.0;

    // CRITICAL FIX: Compute baseline bathymetry from age (missing from unified pipeline!)
    // This was causing geological processes to overwrite each other instead of being additive
    let n = world.grid.cells;
    let previous_depth = world.depth_m.clone(); // Store for structural volume change calculation
    
    // Recompute baseline bathymetry from age (like world.rs does)
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
    
    // Now prepare additive tectonic edits on top of this baseline
    let base_depth = world.depth_m.clone();
    let th_before = world.th_c_m.clone();
    // Prepare persistent scratch buffers
    if world.scratch.f32_a.len() != n {
        world.scratch.f32_a.resize(n, 0.0);
    }
    if world.scratch.f32_b.len() != n {
        world.scratch.f32_b.resize(n, 0.0);
    }
    if world.scratch.f32_c.len() != n {
        world.scratch.f32_c.resize(n, 0.0);
    }
    if world.scratch.f32_d.len() != n {
        world.scratch.f32_d.resize(n, 0.0);
    }
    let delta_tect: &mut [f32] = &mut world.scratch.f32_a;
    for v in delta_tect.iter_mut() {
        *v = 0.0;
    }
    let sub_delta_step: &mut [f32] = &mut world.scratch.f32_b; // hold CFL-limited subduction deltas
    for v in sub_delta_step.iter_mut() {
        *v = 0.0;
    }
    let tr_delta_step: &mut [f32] = &mut world.scratch.f32_c; // hold CFL-limited transform deltas
    for v in tr_delta_step.iter_mut() {
        *v = 0.0;
    }
    let cap_per_step: f32 = 200.0; // global trust region cap (meters)
                                   // Consolidated CFL summaries and staged deltas for exclusivity resolution
    let mut cfl_summaries: Vec<String> = Vec::new();
    // Reuse scratch for CFL-limited deltas
    let sub_delta_step: &mut [f32] = &mut world.scratch.f32_b;
    for v in sub_delta_step.iter_mut() {
        *v = 0.0;
    }
    let tr_delta_step: &mut [f32] = &mut world.scratch.f32_c;
    for v in tr_delta_step.iter_mut() {
        *v = 0.0;
    }
    // DL-1 metrics accumulators; assigned deterministically before use
    let dl1_composite_caps: u32;
    let dl1_th_caps: u32;
    let dl1_rift_soft_cells: u32;
    let dl1_rift_uplift_soft_cells: u32;
    // DL-2 mass/volume budget accumulators (m^3)
    let mut dl2_dv_sub_m3: f64 = 0.0;
    if cfg.enable_subduction && do_sub {
        let t_s0 = std::time::Instant::now();
        let _span = tracing::info_span!("subduction").entered();
        // Unified CFL limiter for subduction band deltas (no automatic substep escalation)
        let sub_delta = &mut world.scratch.f32_d; // reuse scratch
        for v in sub_delta.iter_mut() {
            *v = 0.0;
        }
        let _sub_stats = subduction::compute_subduction_delta(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            sub_delta,
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

        // Characteristic width: trench half-width (meters)
        let width_m = (cfg.sub_trench_half_width_km as f64).max(1.0) * 1000.0;
        let mut cfl_stats = crate::cfl::CflStats::default();
        let cfl_cfg = crate::cfl::CflConfig { max_cfl: 0.3, debug_log: false };
        {
            let threads =
                thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let sub_slice = sub_delta[start..end].to_vec();
                    let area = world.area_m2[start..end].to_vec();
                    let cfl_cfg_local = cfl_cfg;
                    let width_m_local = width_m;
                    let dt_local = dt as f64;
                    let handle = scope.spawn(move || {
                        let mut dv_local: f64 = 0.0;
                        let mut s = crate::cfl::CflStats::default();
                        let mut out_local = vec![0.0f32; sub_slice.len()];
                        for k in 0..sub_slice.len() {
                            let res = crate::cfl::limit(
                                sub_slice[k] as f64,
                                width_m_local,
                                dt_local,
                                cfl_cfg_local,
                            );
                            s.update(&res);
                            out_local[k] = res.scaled_displacement_m as f32;
                            if res.scaled_displacement_m > 0.0 {
                                dv_local += (res.scaled_displacement_m as f64) * (area[k] as f64);
                            }
                        }
                        (start, s, out_local, dv_local)
                    });
                    let (start, s, out_local, dv_local) = handle.join().unwrap();
                    let len = out_local.len();
                    sub_delta_step[start..start + len].copy_from_slice(&out_local);
                    dl2_dv_sub_m3 += dv_local;
                    cfl_stats.cells_processed += s.cells_processed;
                    cfl_stats.cells_scaled += s.cells_scaled;
                    if s.max_cfl > cfl_stats.max_cfl {
                        cfl_stats.max_cfl = s.max_cfl;
                    }
                    if s.min_scale_factor < cfl_stats.min_scale_factor {
                        cfl_stats.min_scale_factor = s.min_scale_factor;
                    }
                    // accumulate mean via weighted sum
                    if cfl_stats.cells_processed > 0 {
                        let prev = cfl_stats.mean_cfl
                            * ((cfl_stats.cells_processed - s.cells_processed) as f64);
                        cfl_stats.mean_cfl = (prev + s.mean_cfl * (s.cells_processed as f64))
                            / (cfl_stats.cells_processed as f64);
                    }
                }
            });
        }
        cfl_summaries.push(format!(
            "subduction max={:.3} mean={:.3} scaled%={:.1}",
            cfl_stats.max_cfl,
            cfl_stats.mean_cfl,
            if cfl_stats.cells_processed > 0 {
                100.0 * (cfl_stats.cells_scaled as f64) / (cfl_stats.cells_processed as f64)
            } else {
                0.0
            }
        ));
        ms_subduction += t_s0.elapsed().as_secs_f64() * 1000.0;
    }
    // Continental rifting (delta-only) before transforms
    {
        let p_rift = crate::rifting::RiftingParams {
            c_rift_min: 0.60,
            v_open_min_m_per_yr: 0.020, // 20 mm/yr
            w_core_km: 60.0,
            w_taper_km: 200.0,
            k_thin: 0.0025,
            alpha_subs: 0.33,
            ocean_thresh: 0.15,
            k_c_oceanize: 0.030,
            reset_age_on_core: true,
            enable_shoulder: true,
            w_bulge_km: 80.0,
            beta_shoulder: 0.10,
            couple_flexure: false,
            thc_min_m: 25_000.0,
            thc_max_m: 65_000.0,
        };
        let rift_tmp = &mut world.scratch.f32_d; // reuse scratch
        for v in rift_tmp.iter_mut() {
            *v = 0.0;
        }
        let r_stats = crate::rifting::apply_rifting(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.plates.kind,
            &mut world.c,
            &mut world.th_c_m,
            &mut world.age_myr,
            rift_tmp,
            &world.area_m2,
            &p_rift,
            dt as f64,
        );
        dl1_rift_soft_cells = r_stats.cells_rate_softened as u32;
        dl1_rift_uplift_soft_cells = r_stats.cells_uplift_softened as u32;
        let mut units_bug_max = 0.0f32;
        {
            let threads =
                thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let src = rift_tmp[start..end].to_vec();
                    let handle = scope.spawn(move || {
                        let mut local_units = 0.0f32;
                        let mut adds = vec![0.0f32; src.len()];
                        for k in 0..src.len() {
                            let raw = src[k];
                            if raw.abs() > 1000.0 {
                                local_units = local_units.max(raw.abs());
                            }
                            adds[k] = raw.clamp(-cap_per_step, cap_per_step);
                        }
                        (start, local_units, adds)
                    });
                    let (s0, lu, adds) = handle.join().unwrap();
                    if lu > units_bug_max {
                        units_bug_max = lu;
                    }
                    for (k, &a) in adds.iter().enumerate() {
                        delta_tect[s0 + k] += a;
                    }
                }
            });
        }
        if units_bug_max > 0.0 {
            println!("[UNITS_BUG] rifting Δz max {:.1} m (>1000 m)", units_bug_max);
        }
    }
    // Transforms after rifting
    if do_trf {
        let t_t0 = std::time::Instant::now();
        let _span = tracing::info_span!("transforms").entered();
        // Use unified CFL limiter instead of substeps
        let mut cfl_stats = crate::cfl::CflStats::default();
        let cfl_config = crate::cfl::CflConfig { max_cfl: 0.3, debug_log: false };

        // RL-3: Characteristic width for transforms W_t (km) → meters.
        // Use a conservative base 50 km clamped to [20, 80] km.
        let width_t_m: f64 = (50.0f64).clamp(20.0, 80.0) * 1000.0;

        let tr_tmp = &mut world.scratch.f32_d; // reuse scratch
        for v in tr_tmp.iter_mut() {
            *v = 0.0;
        }
        let _ = transforms::apply_transforms(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.v_en,
            tr_tmp,
            transforms::TransformParams {
                tau_open_m_per_yr: 0.005,
                min_tangential_m_per_yr: 0.015,
                max_normal_m_per_yr: 0.001,
                basin_half_width_km: 10.0,
                ridge_like_uplift_m: -3.0,
                basin_deepen_m: 6.0,
            },
            dt as f64,
        );

        // Apply CFL limiting to each cell's displacement
        {
            let threads =
                thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            let mut total_proc = 0u32;
            let mut total_scaled = 0u32;
            let mut max_cfl_all = 0.0f64;
            let mut min_scale_all = 1.0f64;
            let mut mean_weighted_sum = 0.0f64;
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let src = tr_tmp[start..end].to_vec();
                    let cfg_local = cfl_config;
                    let width_local = width_t_m;
                    let dt_local = dt as f64;
                    let handle = scope.spawn(move || {
                        let mut s = crate::cfl::CflStats::default();
                        let mut out_local = vec![0.0f32; src.len()];
                        for k in 0..src.len() {
                            let res =
                                crate::cfl::limit(src[k] as f64, width_local, dt_local, cfg_local);
                            s.update(&res);
                            out_local[k] = res.scaled_displacement_m as f32;
                        }
                        (start, s, out_local)
                    });
                    let (s0, s, out_local) = handle.join().unwrap();
                    let len = out_local.len();
                    tr_delta_step[s0..s0 + len].copy_from_slice(&out_local);
                    total_proc += s.cells_processed;
                    total_scaled += s.cells_scaled;
                    if s.max_cfl > max_cfl_all {
                        max_cfl_all = s.max_cfl;
                    }
                    if s.min_scale_factor < min_scale_all {
                        min_scale_all = s.min_scale_factor;
                    }
                    mean_weighted_sum += (s.mean_cfl as f64) * (s.cells_processed as f64);
                }
            });
            if total_proc > 0 {
                cfl_stats.cells_processed = total_proc;
                cfl_stats.cells_scaled = total_scaled;
                cfl_stats.max_cfl = max_cfl_all;
                cfl_stats.min_scale_factor = min_scale_all;
                cfl_stats.mean_cfl = mean_weighted_sum / (total_proc as f64);
            }
        }

        // RL-4: accumulate consolidated summary instead of per-process prints
        cfl_summaries.push(format!(
            "transforms max={:.3} mean={:.3} scaled%={:.1}",
            cfl_stats.max_cfl,
            cfl_stats.mean_cfl,
            if cfl_stats.cells_processed > 0 {
                100.0 * (cfl_stats.cells_scaled as f64) / (cfl_stats.cells_processed as f64)
            } else {
                0.0
            }
        ));
        ms_transforms += t_t0.elapsed().as_secs_f64() * 1000.0;
    }
    // Continents buoyancy: compute target amplitude directly each step and compose delta relative to last amplitude
    {
        let t_bu0 = std::time::Instant::now();
        // Densities
        let pc = crate::PhysConsts::default();
        let rho_m = pc.rho_m_kg_per_m3;
        let rho_c = pc.rho_c_kg_per_m3;
        let rho_w = pc.rho_w_kg_per_m3;
        let rho_a = pc.rho_air_kg_per_m3;
        let mut dstep_max = 0.0f32;
        let mut amp_max = 0.0f32;
        let mut caps_world_nonzero = 0u32;
        let _step_cap = cap_per_step; // reuse global per-step cap (retained for context)
        if world.delta_buoy_m.len() != n {
            world.delta_buoy_m.resize(n, 0.0);
        }
        // Count non-zero caps (thickness above small epsilon)
        for i in 0..n {
            if world.th_c_m[i] > 1.0 && world.c[i] > 0.0 {
                caps_world_nonzero += 1;
            }
        }
        let mut clamp_count: u32 = 0;
        for (i, dt) in delta_tect.iter_mut().enumerate().take(n) {
            let c = world.c[i].clamp(0.0, 1.0);
            let th = world.th_c_m[i].clamp(25_000.0, 65_000.0);
            let elev_prev = world.sea.eta_m - world.depth_m[i];
            let rho_top = if elev_prev > 0.0 { rho_a } else { rho_w };
            let frac = (rho_m - rho_c) / (rho_m - rho_top);
            let target = -c * frac * (th - 35_000.0);
            let old = world.delta_buoy_m[i];
            let raw_step = target - old;
            // Trust-region clamp consistent with other per-step edits
            let step = raw_step.clamp(-cap_per_step, cap_per_step);
            if step != raw_step {
                clamp_count += 1;
            }
            world.delta_buoy_m[i] = old + step;
            dstep_max = dstep_max.max(step.abs());
            amp_max = amp_max.max(target.abs());
            *dt += step;
        }
        println!("[caps] world_nonzero={}", caps_world_nonzero);
        if clamp_count > 0 {
            let pct = (100.0 * (clamp_count as f64) / (n as f64)) as f32;
            println!(
                "[buoyancy] dZ_step max={:.1} m | amp max={:.0} m | clamped={} ({:.1}%)",
                dstep_max, amp_max, clamp_count, pct
            );
        } else {
            println!("[buoyancy] dZ_step max={:.1} m | amp max={:.0} m", dstep_max, amp_max);
        }
        world.epoch_continents = world.epoch_continents.wrapping_add(1);
        ms_buoyancy += t_bu0.elapsed().as_secs_f64() * 1000.0;
    }
    // Accretion and Orogeny deltas after continents uplift
    {
        let t_ao0 = std::time::Instant::now();
        // Reuse last subduction masks by recomputing quickly if needed
        let sub_tmp = &mut world.scratch.f32_d; // reuse scratch
        for v in sub_tmp.iter_mut() {
            *v = 0.0;
        }
        let sub = subduction::compute_subduction_delta(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            sub_tmp,
            subduction::SubductionParams {
                tau_conv_m_per_yr: (cfg.sub_tau_conv_m_per_yr.max(0.020)) as f64,
                trench_half_width_km: cfg.sub_trench_half_width_km as f64,
                arc_offset_km: cfg.sub_arc_offset_km as f64,
                arc_half_width_km: cfg.sub_arc_half_width_km as f64,
                backarc_width_km: cfg.sub_backarc_width_km as f64,
                trench_deepen_m: cfg.sub_trench_deepen_m,
                arc_uplift_m: cfg.sub_arc_uplift_m.min(8.0),
                backarc_uplift_m: cfg.sub_backarc_uplift_m.min(6.0),
                rollback_offset_m: cfg.sub_rollback_offset_m as f64,
                rollback_rate_km_per_myr: cfg.sub_rollback_rate_km_per_myr as f64,
                backarc_extension_mode: cfg.sub_backarc_extension_mode,
                backarc_extension_deepen_m: cfg.sub_backarc_extension_deepen_m,
                continent_c_min: cfg.sub_continent_c_min,
            },
            Some(&world.c),
        );
        // Exclusive masks: build transform masks and strip out ridge/subduction
        let tr_raw = &mut world.scratch.f32_d; // reuse scratch
        for v in tr_raw.iter_mut() {
            *v = 0.0;
        }
        let (tr_masks, _tr_stats) = transforms::apply_transforms(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.v_en,
            tr_raw,
            transforms::TransformParams {
                tau_open_m_per_yr: 0.005,
                min_tangential_m_per_yr: 0.015,
                max_normal_m_per_yr: 0.001,
                basin_half_width_km: 10.0,
                ridge_like_uplift_m: -3.0,
                basin_deepen_m: 6.0,
            },
            dt as f64,
        );
        // Build exclusive masks (binary); simple 1-cell dilation via neighbor OR
        let mut m_s: Vec<u8> = vec![0; n];
        for (i, m) in m_s.iter_mut().enumerate().take(n) {
            if sub.masks.trench[i] || sub.masks.arc[i] || sub.masks.backarc[i] {
                *m = 1;
            }
        }
        let mut m_t: Vec<u8> = vec![0; n];
        for (i, m) in m_t.iter_mut().enumerate().take(n) {
            if tr_masks.pull_apart[i] || tr_masks.restraining[i] {
                *m = 1;
            }
        }
        let m_t_before = m_t.clone();
        // Dilate subduction mask by 1 ring
        let mut m_s_dil = m_s.clone();
        for (i, &ms) in m_s.iter().enumerate().take(n) {
            if ms != 0 {
                for &nj in &world.grid.n1[i] {
                    m_s_dil[nj as usize] = 1;
                }
            }
        }
        // Exclude ridge (not explicitly built here) by leaving room for future ridge mask
        // Exclusivity: strip subduction from transforms
        for i in 0..n {
            if m_s_dil[i] != 0 {
                m_t[i] = 0;
            }
        }
        // Diagnostics: overlaps before resolution
        let mut overlap_before = 0u32;
        for i in 0..n {
            if m_s_dil[i] != 0 && m_t_before[i] != 0 {
                overlap_before += 1;
            }
        }
        if overlap_before > 0 {
            println!(
                "[excl] subduction∩transform candidates before resolve: {} cells",
                overlap_before
            );
        }
        // Apply deltas honoring exclusivity
        let mut overlaps_applied = 0u32;
        {
            let threads =
                std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let ms = m_s_dil[start..end].to_vec();
                    let mt = m_t[start..end].to_vec();
                    let sub = sub_delta_step[start..end].to_vec();
                    let trd = tr_delta_step[start..end].to_vec();
                    let handle = scope.spawn(move || {
                        let mut local_overlap = 0u32;
                        let mut adds = vec![0.0f32; ms.len()];
                        for k in 0..ms.len() {
                            let mut applied = 0u8;
                            if ms[k] != 0 {
                                adds[k] += sub[k].clamp(-cap_per_step, cap_per_step);
                                applied |= 1;
                            }
                            if mt[k] != 0 {
                                adds[k] += trd[k].clamp(-cap_per_step, cap_per_step);
                                applied |= 2;
                            }
                            if applied == 3 {
                                local_overlap += 1;
                            }
                        }
                        (start, local_overlap, adds)
                    });
                    let (s0, ov, adds) = handle.join().unwrap();
                    overlaps_applied += ov;
                    for (k, &a) in adds.iter().enumerate() {
                        delta_tect[s0 + k] += a;
                    }
                }
            });
        }
        if overlaps_applied > 0 {
            println!("[excl][ERROR] stacked edits applied to {} cells", overlaps_applied);
        }
        debug_assert_eq!(overlaps_applied, 0, "Exclusivity failed: stacked edits applied");
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
        let acc_tmp = &mut world.scratch.f32_d; // reuse scratch
        for v in acc_tmp.iter_mut() {
            *v = 0.0;
        }
        let _ = crate::accretion::apply_oc_accretion(
            &world.grid,
            &sub.masks,
            &world.boundaries,
            &vel3,
            &world.plates.kind,
            &mut world.c,
            &mut world.th_c_m,
            acc_tmp,
            &world.area_m2,
            &p_acc,
            dt as f64,
        );
        let mut units_bug_max_acc = 0.0f32;
        {
            let threads =
                std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let src = acc_tmp[start..end].to_vec();
                    let handle = scope.spawn(move || {
                        let mut local_units = 0.0f32;
                        let mut adds = vec![0.0f32; src.len()];
                        for k in 0..src.len() {
                            let raw = src[k];
                            if raw.abs() > 1000.0 {
                                local_units = local_units.max(raw.abs());
                            }
                            adds[k] = raw.clamp(-cap_per_step, cap_per_step);
                        }
                        (start, local_units, adds)
                    });
                    let (s0, lu, adds) = handle.join().unwrap();
                    if lu > units_bug_max_acc {
                        units_bug_max_acc = lu;
                    }
                    for (k, &a) in adds.iter().enumerate() {
                        delta_tect[s0 + k] += a;
                    }
                }
            });
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
        let orog_tmp = &mut world.scratch.f32_d; // reuse scratch
        for v in orog_tmp.iter_mut() {
            *v = 0.0;
        }
        let _ = crate::orogeny::apply_cc_orogeny(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &vel3,
            &world.c,
            &world.area_m2,
            &mut world.th_c_m,
            orog_tmp,
            &p_orog,
            dt as f64,
            &world.plates.kind,
        );
        let mut units_bug_max_orog = 0.0f32;
        {
            let threads =
                std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n + threads - 1) / threads).max(1);
            std::thread::scope(|scope| {
                for t in 0..threads {
                    let start = t * chunk;
                    if start >= n {
                        break;
                    }
                    let end = (start + chunk).min(n);
                    let src = orog_tmp[start..end].to_vec();
                    let handle = scope.spawn(move || {
                        let mut local_units = 0.0f32;
                        let mut adds = vec![0.0f32; src.len()];
                        for k in 0..src.len() {
                            let raw = src[k];
                            if raw.abs() > 1000.0 {
                                local_units = local_units.max(raw.abs());
                            }
                            adds[k] = raw.clamp(-300.0, 300.0);
                        }
                        (start, local_units, adds)
                    });
                    let (s0, lu, adds) = handle.join().unwrap();
                    if lu > units_bug_max_orog {
                        units_bug_max_orog = lu;
                    }
                    for (k, &a) in adds.iter().enumerate() {
                        delta_tect[s0 + k] += a;
                    }
                }
            });
        }
        if units_bug_max_orog > 0.0 {
            println!("[UNITS_BUG] orogeny Δz max {:.1} m (>1000 m)", units_bug_max_orog);
        }

        // Composite per-cell cap on the sum across operators
        let mut composite_caps = 0u32;
        // Single-thread clamp to avoid cross-thread mutation complexity
        for v in delta_tect.iter_mut() {
            if v.abs() > cap_per_step {
                *v = v.clamp(-cap_per_step, cap_per_step);
                composite_caps += 1;
            }
        }
        if composite_caps > 0 {
            println!(
                "[cap] composite per-cell: clamped {} cells (cap={:.1} m)",
                composite_caps, cap_per_step
            );
        }
        dl1_composite_caps = composite_caps;

        // Global trust region: scale deltas if any exceed per-step hard limit
        let mut max_abs = 0.0f32;
        for &d in delta_tect.iter() {
            max_abs = max_abs.max(d.abs());
        }
        if max_abs > cap_per_step {
            let scale = cap_per_step / max_abs;
            for d in delta_tect.iter_mut() {
                *d *= scale;
            }
            println!("[trust] scaled tectonic deltas by {:.3} (max {:.1} m)", scale, max_abs);
        }
        // Enforce per-step cap on crustal thickness changes
        let mut th_caps = 0u32;
        let mut th_units_bug = 0.0f32;
        {
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
        }
        if th_caps > 0 {
            println!("[cap] th_c per-step: clamped {} cells (cap={:.1} m)", th_caps, cap_per_step);
        }
        dl1_th_caps = th_caps;
        if th_units_bug > 0.0 {
            println!("[UNITS_BUG] Δth_c max {:.1} m (>1000 m)", th_units_bug);
        }
        ms_accretion_orogeny += t_ao0.elapsed().as_secs_f64() * 1000.0;
    }
    // RL-4: emit one consolidated CFL line per step
    if !cfl_summaries.is_empty() {
        println!("[CFL] {}", cfl_summaries.join(" | "));
    }
    // Build age-only baselines now (after deltas) to free scratch earlier, then compose
    let n_cells = world.grid.cells;
    let dt_f = dt as f64;
    // Reuse scratch for d0_prev
    let d0_prev = &mut world.scratch.f32_d;
    for v in d0_prev.iter_mut() {
        *v = 0.0;
    }
    let mut d0_now = vec![0.0f32; n_cells];
    {
        let threads = thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
        let chunk = ((n_cells + threads - 1) / threads).max(1);
        let mut handles = Vec::new();
        for t in 0..threads {
            let start = t * chunk;
            if start >= n_cells {
                break;
            }
            let end = (start + chunk).min(n_cells);
            // Copy read-only slices for this chunk into owned buffers for the worker
            let ages_now: Vec<f32> = world.age_myr[start..end].to_vec();
            let ages_reset: Vec<f32> = ages[start..end].to_vec();
            handles.push(thread::spawn(move || {
                let len = ages_now.len();
                let mut prev_local = vec![0.0f32; len];
                let mut now_local = vec![0.0f32; len];
                for k in 0..len {
                    let age_now = ages_now[k] as f64;
                    let age_was =
                        if ages_reset[k] == 0.0 { 0.0 } else { (age_now - dt_f).max(0.0) };
                    prev_local[k] =
                        crate::age::depth_from_age_blend(age_was, 2500.0, 350.0, 6300.0, 60.0)
                            as f32;
                    now_local[k] =
                        crate::age::depth_from_age_blend(age_now, 2500.0, 350.0, 6300.0, 60.0)
                            as f32;
                }
                (start, prev_local, now_local)
            }));
        }
        for h in handles {
            if let Ok((start, prev_local, now_local)) = h.join() {
                let len = prev_local.len();
                d0_prev[start..start + len].copy_from_slice(&prev_local);
                d0_now[start..start + len].copy_from_slice(&now_local);
            }
        }
    }
    // Preserve previous structural component relative to age baseline, then add incremental deltas
    {
        let threads = thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
        let chunk = ((n + threads - 1) / threads).max(1);
        // Use owned chunks to avoid borrowing `world.depth_m` across scoped threads
        let mut out_chunks: Vec<(usize, Vec<f32>)> = Vec::new();
        std::thread::scope(|scope| {
            for t in 0..threads {
                let start = t * chunk;
                if start >= n {
                    break;
                }
                let end = (start + chunk).min(n);
                let base = base_depth[start..end].to_vec();
                let dprev = d0_prev[start..end].to_vec();
                let dnow = d0_now[start..end].to_vec();
                let dtect = delta_tect[start..end].to_vec();
                let handle = scope.spawn(move || {
                    let mut out_local = vec![0.0f32; base.len()];
                    for k in 0..out_local.len() {
                        let struct_prev = base[k] - dprev[k];
                        out_local[k] = (dnow[k] + struct_prev + dtect[k]).clamp(-8000.0, 8000.0);
                    }
                    (start, out_local)
                });
                let (s0, out_local) = handle.join().unwrap();
                out_chunks.push((s0, out_local));
            }
        });
        for (s0, out_local) in out_chunks {
            let len = out_local.len();
            world.depth_m[s0..s0 + len].copy_from_slice(&out_local);
        }
    }

    // Flexure: experimental GPU V-cycle or Winkler fallback
    if cfg.enable_flexure && do_flx {
        let t_fx0 = std::time::Instant::now();
        // Update Te field from age/C/th_c before assembling load
        crate::flexure::compute_te_field(&world.age_myr, &world.c, &world.th_c_m, &mut world.te_m);
        let pc = crate::PhysConsts::default();
        let lp = flexure_loads::LoadParams {
            rho_w: pc.rho_w_kg_per_m3,
            rho_c: pc.rho_c_kg_per_m3,
            g: pc.g_m_per_s2,
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
        // Depth sanity after flexure
        {
            let mut bad = 0u32;
            let mut dmin = f32::INFINITY;
            let mut dmax = f32::NEG_INFINITY;
            let mut sum = 0.0f64;
            let mut n = 0u32;
            for &d in world.depth_m.iter() {
                if !d.is_finite() {
                    bad += 1;
                    continue;
                }
                n += 1;
                if d < dmin {
                    dmin = d;
                }
                if d > dmax {
                    dmax = d;
                }
                sum += d as f64;
            }
            if bad > 0
                || !dmin.is_finite()
                || !dmax.is_finite()
                || dmin >= dmax
                || dmin.abs() > 20_000.0
                || dmax.abs() > 20_000.0
            {
                println!(
                    "[sanity][after flexure] n_bad={} min/mean/max={:.0}/{:.0}/{:.0} m",
                    bad,
                    dmin,
                    if n > 0 { (sum / (n as f64)) as f32 } else { 0.0 },
                    dmax
                );
            }
        }
        ms_flexure += t_fx0.elapsed().as_secs_f64() * 1000.0;
    }
    // NOTE(Science): Winkler foundation ignores lithospheric flexural coupling length scales.
    // It is useful as a placeholder but cannot reproduce forebulge/trough patterns.

    // Erosion/diffusion
    let do_surf = k % (cfg.cadence_surf_every.max(1) as u64) == 0;
    if cfg.enable_erosion && do_surf {
        let t_sf0 = std::time::Instant::now();
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
        // Depth sanity after surface processes
        {
            let mut bad = 0u32;
            let mut dmin = f32::INFINITY;
            let mut dmax = f32::NEG_INFINITY;
            let mut sum = 0.0f64;
            let mut n = 0u32;
            for &d in world.depth_m.iter() {
                if !d.is_finite() {
                    bad += 1;
                    continue;
                }
                n += 1;
                if d < dmin {
                    dmin = d;
                }
                if d > dmax {
                    dmax = d;
                }
                sum += d as f64;
            }
            if bad > 0
                || !dmin.is_finite()
                || !dmax.is_finite()
                || dmin >= dmax
                || dmin.abs() > 20_000.0
                || dmax.abs() > 20_000.0
            {
                println!(
                    "[sanity][after surface] n_bad={} min/mean/max={:.0}/{:.0}/{:.0} m",
                    bad,
                    dmin,
                    if n > 0 { (sum / (n as f64)) as f32 } else { 0.0 },
                    dmax
                );
            }
        }
        ms_surface += t_sf0.elapsed().as_secs_f64() * 1000.0;
    }

    // Depth sanity guard: reject non-finite or absurd magnitudes; revert to previous depth
    {
        let n = world.depth_m.len().min(base_depth.len());
        for (i, _) in base_depth.iter().enumerate().take(n) {
            let d = world.depth_m[i];
            if !d.is_finite() || d.abs() > 20_000.0 {
                world.depth_m[i] = base_depth[i];
            }
        }
    }

    // CRITICAL FIX: Apply continental uplift AFTER all geological processes finish
    // This ensures continental elevation follows plate motion and isn't overwritten
    if cfg.enable_rigid_motion {
        let t_uplift0 = std::time::Instant::now();
        crate::continent::apply_uplift_from_c_thc(&mut world.depth_m, &world.c, &world.th_c_m);
        let ms_continental_uplift = t_uplift0.elapsed().as_secs_f64() * 1000.0;
        println!("[continental_uplift] Applied continental thickness to depth_m field AFTER geological processes ({:.2} ms)", ms_continental_uplift);
    }

    // Solve eta to hit target land fraction unless frozen; elevation is independent of eta
    // Eta solve moved after writing elevation to reduce coastline speckle

    if !cfg.freeze_eta && do_sea {
        let t_eta0 = std::time::Instant::now();
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
        if let Some(_r) = world.sea_level_ref {
            // Canonical land fraction/eta solve on CPU (area-weighted) with guard rails
            // Sanitize/scan inputs
            // scan inputs; avoid explicit counters for clippy cleanliness
            let mut n_bad = 0u32;
            let mut dmin = f32::INFINITY;
            let mut dmax = f32::NEG_INFINITY;
            for &d in world.depth_m.iter() {
                if d.is_finite() {
                    if d < dmin {
                        dmin = d;
                    }
                    if d > dmax {
                        dmax = d;
                    }
                } else {
                    n_bad += 1;
                }
            }
            // If inputs are suspect, hold previous eta
            if n_bad == 0 && dmin.is_finite() && dmax.is_finite() && dmin < dmax {
                let eta_prev = *surf.eta_m as f64;
                let eta_candidate = crate::isostasy::solve_eta_for_land(
                    &world.depth_m,
                    &world.area_m2,
                    cfg.target_land_frac as f64,
                    1e-4,
                    64,
                );
                // Validate and clamp per-frame change
                let soft_lo = (dmin as f64) - 10_000.0;
                let soft_hi = (dmax as f64) + 10_000.0;
                if eta_candidate.is_finite() && eta_candidate >= soft_lo && eta_candidate <= soft_hi
                {
                    // Pre-clamp diagnostic at candidate eta (tests solver itself)
                    let land_at_candidate = crate::isostasy::land_fraction_with_eta(
                        &world.depth_m,
                        &world.area_m2,
                        eta_candidate as f32,
                    );
                    let topo_frozen = !(cfg.enable_rigid_motion
                        || cfg.enable_subduction
                        || cfg.enable_erosion
                        || cfg.enable_flexure);
                    if topo_frozen {
                        debug_assert!(
                            (land_at_candidate as f32 - cfg.target_land_frac).abs() <= 0.002 + 1e-6,
                            "land mismatch at candidate: land={:.3} target={:.3}",
                            land_at_candidate,
                            cfg.target_land_frac
                        );
                    }
                    let dmax_step = 200.0_f64 * (dt as f64); // m per frame (scaled by dt Myr)
                    let d_eta = (eta_candidate - eta_prev).clamp(-dmax_step, dmax_step);
                    let eta_final = (eta_prev + d_eta).clamp(soft_lo, soft_hi);
                    *surf.eta_m = eta_final as f32;
                    // Diagnostics: report both candidate and final land fractions
                    let land_final = crate::isostasy::land_fraction_with_eta(
                        &world.depth_m,
                        &world.area_m2,
                        *surf.eta_m,
                    );
                    println!(
                        "[sea] target_land={:.3} eta_cand={:+.0} m land_cand={:.3} | eta={:+.0} m land_cpu={:.3}",
                        cfg.target_land_frac as f64,
                        eta_candidate,
                        land_at_candidate,
                        *surf.eta_m as f64,
                        land_final
                    );
                    // Assert close to target only when topology is frozen
                    if topo_frozen {
                        debug_assert!(
                            (land_final as f32 - cfg.target_land_frac).abs() <= 0.002 + 1e-6,
                            "land mismatch: land_cpu={:.3} target={:.3}",
                            land_final,
                            cfg.target_land_frac
                        );
                    }
                } else {
                    // Hold previous eta
                    // Optionally log once per second, here we keep silent to avoid spam
                }
            } else {
                // Hold previous eta if depths had non-finites
            }
        }
        ms_eta += t_eta0.elapsed().as_secs_f64() * 1000.0;
    }

    // Write surface elevation after eta solve: z = eta - depth
    for (e, &d) in surf.elevation_m.iter_mut().zip(world.depth_m.iter()) {
        let z = (*surf.eta_m as f64 - d as f64) as f32;
        *e = z;
    }

    // Invariants & guards (TST-2): check after we have elevation
    check_invariants(world, Some(surf.elevation_m));

    // Optional: mass/volume budgets and per-step deltas (DL-2)
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
        // DL-2 additional budget terms
        let mut dv_th_m3: f64 = 0.0;
        let mut dv_th_pos_m3: f64 = 0.0;
        let mut dv_th_neg_m3: f64 = 0.0;
        for (i, &a_cell) in world.area_m2.iter().enumerate().take(n) {
            let dth = (world.th_c_m[i] - th_before[i]) as f64;
            let dv = (a_cell as f64) * dth;
            dv_th_m3 += dv;
            if dv >= 0.0 {
                dv_th_pos_m3 += dv;
            } else {
                dv_th_neg_m3 += dv;
            }
        }
        let mut dv_struct_m3: f64 = 0.0;
        for (i, &a) in world.area_m2.iter().enumerate().take(n) {
            let dd = (world.depth_m[i] as f64) - (previous_depth[i] as f64);
            dv_struct_m3 += (a as f64) * dd;
        }
        let pc = crate::PhysConsts::default();
        let k_b =
            (pc.rho_m_kg_per_m3 as f64 - pc.rho_c_kg_per_m3 as f64) / (pc.rho_m_kg_per_m3 as f64);
        let dv_buoy_m3 = k_b * dv_th_m3;
        let closure = dv_struct_m3 - dv_buoy_m3;
        let tol = 0.05f64;
        let closure_ok = if dv_buoy_m3.abs() > 0.0 {
            (closure.abs() / dv_buoy_m3.abs()) <= tol
        } else {
            closure.abs() <= 1.0
        };
        println!(
            "[budget] land={:.1}% M_cont={:.3e} (d={:+.1e}) M_sed={:.3e} (d={:+.1e}) V_ocean={:.3e} (d={:+.1e}) err_ref={:+.1e} | dV_th={:+.3e} (pos={:+.3e} neg={:+.3e}) dV_sub={:+.3e} dV_struct={:+.3e} dV_buoy={:+.3e} closure={:+.3e} {}",
            100.0 * land_frac,
            m_cont,
            d_cont,
            m_sed,
            d_sed,
            vol_ocean,
            d_vol,
            vol_err,
            dv_th_m3,
            dv_th_pos_m3,
            dv_th_neg_m3,
            dl2_dv_sub_m3,
            dv_struct_m3,
            dv_buoy_m3,
            closure,
            if closure_ok { "[ok]" } else { "[WARN]" }
        );
        world.last_mass_cont_kg = m_cont;
        world.last_mass_sed_kg = m_sed;
        world.last_ocean_vol_m3 = vol_ocean;
    }
    // RL-4: emit one consolidated CFL line per step (if any) (throttled: once per step)
    if !cfl_summaries.is_empty() {
        tracing::info!(target: "CFL", summary = %cfl_summaries.join(" | "));
    }
    // DL-1: boundary composition & cap metrics (once per ~Myr)
    if true {
        // Boundary composition by class (edge count proxy)
        let mut n_r = 0u32;
        let mut n_s = 0u32;
        let mut n_t = 0u32;
        for &(_u, _v, cls) in &world.boundaries.edges {
            match cls {
                1 => n_r += 1,
                2 => n_s += 1,
                3 => n_t += 1,
                _ => {}
            }
        }
        let tot = (n_r + n_s + n_t).max(1) as f64;
        let pr = 100.0 * (n_r as f64) / tot;
        let ps = 100.0 * (n_s as f64) / tot;
        let pt = 100.0 * (n_t as f64) / tot;
        // Oceanized cores proxy
        let mut oceanized = 0u32;
        for &cval in &world.c {
            if cval < 0.15 {
                oceanized += 1;
            }
        }
        tracing::info!(target: "diag.boundary", pr = pr, ps = ps, pt = pt, oceanized = oceanized, cap_composite = dl1_composite_caps, cap_thc = dl1_th_caps, rift_soft_rate = dl1_rift_soft_cells, rift_soft_uplift = dl1_rift_uplift_soft_cells);
    }
    // Advance clock for viewer pipeline path
    world.clock.t_myr += dt as f64;
    world.clock.step_idx = world.clock.step_idx.saturating_add(1);

    // Consolidated perf log for viewer pipeline
    let step_ms = ms_boundaries
        + ms_kinematics
        + ms_force_balance
        + ms_ridge
        + ms_subduction
        + ms_transforms
        + ms_buoyancy
        + ms_accretion_orogeny
        + ms_flexure
        + ms_surface
        + ms_eta;
    println!(
        "[perf.pipeline] step={:.2} ms | kin={:.2} | bounds={:.2} | fb={:.2} | ridge={:.2} | subd={:.2} | trf={:.2} | buoy={:.2} | acc+orog={:.2} | flex={:.2} | surf={:.2} | eta={:.2}",
        step_ms,
        ms_kinematics,
        ms_boundaries,
        ms_force_balance,
        ms_ridge,
        ms_subduction,
        ms_transforms,
        ms_buoyancy,
        ms_accretion_orogeny,
        ms_flexure,
        ms_surface,
        ms_eta
    );
}
