//! Centralized physics pipeline (T-632): run full step sequence and expose a single surface elevation.
//!
//! This viewer-facing pipeline contrasts with `world::step_once` in one key way:
//! - Here we compute and output `eta` (sea level) but do not bake the offset back into `world.depth_m`.
//!   Callers must render elevation as `z = -depth - eta`. This reduces coastline flicker and keeps
//!   `depth_m` as a separable tectonic + age baseline.
//! - Subduction/transform/rifting/orogeny still edit `depth_m` additively on top of age-derived bathy.
//!   If these are applied at high cadence without normalization, long-run drift of `depth_m` can occur.

use crate::{boundaries, flexure_loads, ridge, sea_level, subduction, world::World};

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
    // NOTE: The baseline above overwrites previous `depth_m` fully. Any persistent tectonic edits
    // must be applied after this line each frame, or tracked in separate fields and composed.
    // Apply subduction band edits in-place on top of baseline
    if cfg.enable_subduction {
        let _sub_stats = subduction::apply_subduction(
            &world.grid,
            &world.boundaries,
            &world.plates.plate_id,
            &world.age_myr,
            &world.v_en,
            &mut world.depth_m,
            subduction::SubductionParams {
                tau_conv_m_per_yr: 0.005,
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
    }

    // Flexure: respond to loads (simple Winkler-like response for MVP)
    if cfg.enable_flexure {
        let lp = flexure_loads::LoadParams {
            rho_w: 1030.0,
            rho_c: 2900.0,
            g: 9.81,
            sea_level_m: *surf.eta_m,
        };
        let f_load = flexure_loads::assemble_load_from_depth(&world.grid, &world.depth_m, &lp);
        let k = 3.0e8f32; // foundation stiffness (N/m^3)
        for (d, &q) in world.depth_m.iter_mut().zip(f_load.iter()) {
            let w = q / k;
            *d = (*d + w).clamp(-8000.0, 8000.0);
        }
    }
    // NOTE(Science): Winkler foundation ignores lithospheric flexural coupling length scales.
    // It is useful as a placeholder but cannot reproduce forebulge/trough patterns.

    // Erosion/diffusion
    if cfg.enable_erosion {
        let sp = crate::surface::SurfaceParams {
            k_stream: 1.0e-6,
            m_exp: 0.5,
            n_exp: 1.0,
            k_diff: 1.0e-2,
            k_tr: 0.0,
            p_exp: 1.0,
            q_exp: 1.0,
            rho_sed: 2500.0,
            min_slope: 1.0e-4,
            subcycles: 1,
            couple_flexure: false,
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

    if !cfg.freeze_eta {
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
}
