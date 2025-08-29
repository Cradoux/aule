//! Centralized physics pipeline (T-632): run full step sequence and expose a single surface elevation.
//!
//! This viewer-facing pipeline contrasts with `world::step_once` in one key way:
//! - Here we compute and output `eta` (sea level) but do not bake the offset back into `world.depth_m`.
//!   Callers must render elevation as `z = -depth - eta`. This reduces coastline flicker and keeps
//!   `depth_m` as a separable tectonic + age baseline.
//! - Subduction/transform/rifting/orogeny still edit `depth_m` additively on top of age-derived bathy.
//!   If these are applied at high cadence without normalization, long-run drift of `depth_m` can occur.

use crate::{boundaries, flexure_loads, ridge, sea_level, subduction, transforms, world::World};

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

    // Kinematics:
    // - Refresh per-cell velocities from plates
    // - Optionally advect plate_id via semi-Lagrangian nearest-neighbour when enabled
    // - Classify boundaries using the refreshed fields
    if cfg.enable_rigid_motion {
        world.v_en =
            crate::plates::velocity_en_m_per_yr(&world.grid, &world.plates, &world.plates.plate_id);
        let vel3 = crate::plates::velocity_field_m_per_yr(
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
    }

    // Boundaries classification
    world.boundaries =
        boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);

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

    // Build baseline depth from age curve
    for i in 0..world.grid.cells {
        let d0 = crate::age::depth_from_age_plate(
            world.age_myr[i] as f64,
            2600.0,
            world.clock.t_myr,
            6000.0,
            1.0e-6,
        ) as f32;
        world.depth_m[i] = d0.clamp(0.0, 6000.0);
    }
    // Compose tectonic edits as a single delta added to baseline
    let n = world.grid.cells;
    let base_depth = world.depth_m.clone();
    let mut delta_tect: Vec<f32> = vec![0.0; n];
    if cfg.enable_subduction && do_sub {
        let mut sub_tmp = vec![0.0f32; n];
        let _sub_stats = subduction::apply_subduction(
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
        for i in 0..n {
            delta_tect[i] += sub_tmp[i] - base_depth[i];
        }
    }
    if do_trf {
        let mut tr_tmp = vec![0.0f32; n];
        let _ = transforms::apply_transforms(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.v_en,
            &mut tr_tmp,
            transforms::TransformParams {
                tau_open_m_per_yr: 0.005,
                min_tangential_m_per_yr: 0.003,
                basin_half_width_km: 50.0,
                ridge_like_uplift_m: -200.0,
                basin_deepen_m: 300.0,
            },
        );
        for i in 0..n {
            delta_tect[i] += tr_tmp[i];
        }
    }
    for i in 0..n {
        world.depth_m[i] = (base_depth[i] + delta_tect[i]).clamp(-8000.0, 8000.0);
    }

    // Flexure: experimental GPU V-cycle or Winkler fallback
    if cfg.enable_flexure && do_flx {
        let lp = flexure_loads::LoadParams {
            rho_w: 1030.0,
            rho_c: 2900.0,
            g: 9.81,
            sea_level_m: *surf.eta_m,
        };
        let mut f_load = flexure_loads::assemble_load_from_depth(&world.grid, &world.depth_m, &lp);
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
                // One-tile approximation over full grid (no halos), dx ≈ mean cell scale
                let mut area_mean: f64 = 0.0;
                for a in &world.area_m2 {
                    area_mean += *a as f64;
                }
                let cell_scale_m = (area_mean.max(1.0) / (n as f64)).sqrt() as f32;
                // Map sphere to a strip tile (approx); future: atlas tiling
                let dims = crate::flexure_gpu::TileDims {
                    width: n as u32,
                    height: 1,
                    halo: 2,
                    dx: cell_scale_m,
                };
                let w_tot = dims.width_total() as usize;
                let h_tot = dims.height_total() as usize;
                let n_tot = w_tot * h_tot;
                let mut w_tex = crate::flexure_gpu::GpuTex::new_storage(ctx, n_tot, "flex.w");
                let f_tex = crate::flexure_gpu::GpuTex::new_storage(ctx, n_tot, "flex.f");
                // Stage f into interior row with halos
                let mut f_stage = vec![0.0f32; n_tot];
                let row = dims.halo as usize;
                let start = row * w_tot + (dims.halo as usize);
                f_stage[start..start + n].copy_from_slice(&f_load[..n]);
                ctx.queue.write_buffer(&f_tex.buf, 0, bytemuck::cast_slice(&f_stage));
                let mut scratch = crate::flexure_gpu::FlexScratch::new(ctx, n_tot);
                let e_pa = 70.0e9f32;
                let nu = 0.25f32;
                let te_m = 25_000.0f32;
                let d = e_pa * te_m.powi(3) / (12.0 * (1.0 - nu * nu));
                let p = crate::flexure_gpu::FlexParams {
                    d,
                    k: 0.0,
                    wj_omega: cfg.gpu_wj_omega.clamp(0.4, 0.95),
                    nu1: 1,
                    nu2: 1,
                    levels: cfg.gpu_flex_levels.max(1),
                };
                // Run requested number of V-cycles
                for _ in 0..cfg.gpu_flex_cycles.max(1) {
                    let _stats = flex.v_cycle(ctx, dims, &mut w_tex, &f_tex, &mut scratch, &p);
                }
                // Read back w and extract interior row
                let w_all =
                    crate::flexure_gpu::blocking_read(&ctx.device, &ctx.queue, &w_tex.buf, n_tot);
                let mut w_out = vec![0.0f32; n];
                w_out[..n].copy_from_slice(&w_all[start..start + n]);
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

    // Write surface elevation relative to geoid: z = -depth
    for (e, &d) in surf.elev_m.iter_mut().zip(world.depth_m.iter()) {
        let z = -(d as f64) as f32;
        *e = if z.is_finite() { z } else { 0.0 };
    }

    if !cfg.freeze_eta && do_sea {
        let eta_m =
            sea_level::solve_eta_on_elevation(surf.elev_m, &world.area_m2, cfg.target_land_frac);
        *surf.eta_m = eta_m;
    }

    // Sanitize
    crate::util::sanitize_elevation(surf.elev_m);

    // Optional: mass budget log
    if cfg.log_mass_budget {
        // Proxy mass terms; continent thickness and sediment thickness are explicit
        let mut m_cont: f64 = 0.0;
        let mut m_sed: f64 = 0.0;
        let mut land_area: f64 = 0.0;
        let mut tot_area: f64 = 0.0;
        for i in 0..world.grid.cells {
            let a = world.area_m2[i] as f64;
            m_cont += a * (world.th_c_m.get(i).copied().unwrap_or(0.0) as f64) * 2900.0f64;
            m_sed += a * (world.sediment_m.get(i).copied().unwrap_or(0.0) as f64) * 1800.0f64;
            tot_area += a;
            let elev = -world.depth_m[i];
            if (elev as f64 - *surf.eta_m as f64) > 0.0 {
                land_area += a;
            }
        }
        let land_frac = if tot_area > 0.0 { land_area / tot_area } else { 0.0 };
        println!("[budget] land={:.1}% cont={:.3e} sed={:.3e}", 100.0 * land_frac, m_cont, m_sed);
    }
    // Advance clock for viewer pipeline path
    world.clock.t_myr += dt as f64;
    world.clock.step_idx = world.clock.step_idx.saturating_add(1);
}
