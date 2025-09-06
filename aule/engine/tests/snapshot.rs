use engine::{self, pipeline, world::World};

fn read_last_dl1_record() -> Option<(u32, u32, u32, u32)> {
    let path = std::path::Path::new("out/dl1_metrics.csv");
    if !path.exists() {
        return None;
    }
    let contents = std::fs::read_to_string(path).ok()?;
    let mut last: Option<&str> = None;
    for line in contents.lines() {
        if line.starts_with("t_myr") {
            continue;
        }
        if !line.trim().is_empty() {
            last = Some(line);
        }
    }
    let line = last?;
    let cols: Vec<&str> = line.split(',').collect();
    if cols.len() < 9 {
        return None;
    }
    let cap_comp: u32 = cols[5].parse().ok()?;
    let cap_thc: u32 = cols[6].parse().ok()?;
    let rift_soft_rate: u32 = cols[7].parse().ok()?;
    let rift_soft_uplift: u32 = cols[8].parse().ok()?;
    Some((cap_comp, cap_thc, rift_soft_rate, rift_soft_uplift))
}

#[test]
fn snapshot_tiny_grid_one_step() {
    // Build tiny world (deterministic)
    let f: u32 = 8; // small grid
    let plates = 3u32;
    let seed: u64 = 2025;
    let mut w = World::new(f, plates, seed);
    let mut elev = vec![0.0f32; w.grid.cells];
    let mut eta: f32 = 0.0;
    let surf = pipeline::SurfaceFields { elevation_m: &mut elev, eta_m: &mut eta };
    // Configure to run subduction/transforms each step; enable flexure off for speed
    let cfg = engine::config::PipelineCfg {
        dt_myr: 1.0,
        steps_per_frame: 1,
        enable_flexure: false,
        enable_erosion: false,
        target_land_frac: 0.3,
        freeze_eta: true,
        log_mass_budget: false,
        enable_subduction: true,
        enable_rigid_motion: true,
        cadence_trf_every: 1,
        cadence_sub_every: 1,
        cadence_flx_every: 8,
        cadence_sea_every: 1,
        cadence_surf_every: 8,
        substeps_transforms: 1,
        substeps_subduction: 1,
        use_gpu_flexure: false,
        gpu_flex_levels: 1,
        gpu_flex_cycles: 1,
        gpu_wj_omega: 0.8,
        subtract_mean_load: false,
        sub_tau_conv_m_per_yr: 0.020,
        sub_trench_half_width_km: 40.0,
        sub_arc_offset_km: 140.0,
        sub_arc_half_width_km: 30.0,
        sub_backarc_width_km: 150.0,
        sub_trench_deepen_m: 20.0,
        sub_arc_uplift_m: -12.0,
        sub_backarc_uplift_m: -10.0,
        sub_rollback_offset_m: 0.0,
        sub_rollback_rate_km_per_myr: 0.0,
        sub_backarc_extension_mode: false,
        sub_backarc_extension_deepen_m: 400.0,
        sub_continent_c_min: 0.6,
        surf_k_stream: 1.0e-6,
        surf_m_exp: 0.5,
        surf_n_exp: 1.0,
        surf_k_diff: 0.01,
        surf_k_tr: 0.1,
        surf_p_exp: 0.5,
        surf_q_exp: 1.0,
        surf_rho_sed: 2600.0,
        surf_min_slope: 1.0e-4,
        surf_subcycles: 1,
        surf_couple_flexure: false,
        cadence_spawn_plate_every: 0,
        cadence_retire_plate_every: 0,
        cadence_force_balance_every: 0,
        fb_gain: 0.0,
        fb_damp_per_myr: 0.0,
        fb_k_conv: 1.0,
        fb_k_div: 1.0,
        fb_k_trans: 1.0,
        fb_max_domega: 1.0,
        fb_max_omega: 1.0,
    };
    // One step
    pipeline::step_full(&mut w, surf, cfg);

    // Snapshot metrics
    // Boundary counts by class
    let mut n_r = 0u32;
    let mut n_s = 0u32;
    let mut n_t = 0u32;
    for &(_u, _v, cls) in &w.boundaries.edges {
        match cls {
            1 => n_r += 1,
            2 => n_s += 1,
            3 => n_t += 1,
            _ => {}
        }
    }
    // Land fraction from elevation
    let mut land_area = 0.0f64;
    let mut tot_area = 0.0f64;
    for (a, &z) in w.area_m2.iter().zip(elev.iter()) {
        tot_area += *a as f64;
        if z > 0.0 {
            land_area += *a as f64;
        }
    }
    let land_frac = if tot_area > 0.0 { land_area / tot_area } else { 0.0 };

    // Cap-hit counts from CSV (composite and th_c); optional if file missing
    let caps = read_last_dl1_record();

    // CFL stats: recompute for subduction and transforms like pipeline
    let mut cfl = engine::cfl::CflStats::default();
    // Subduction raw deltas
    let mut sub_tmp = vec![0.0f32; w.grid.cells];
    let _ = engine::subduction::compute_subduction_delta(
        &w.grid,
        &w.boundaries,
        &w.plates.plate_id,
        &w.age_myr,
        &w.v_en,
        &mut sub_tmp,
        engine::subduction::SubductionParams {
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
        Some(&w.c),
    );
    let width_sub_m = (cfg.sub_trench_half_width_km as f64).max(1.0) * 1000.0;
    let cfl_cfg = engine::cfl::CflConfig { max_cfl: 0.3, debug_log: false };
    for &raw in &sub_tmp {
        let res = engine::cfl::limit(raw as f64, width_sub_m, cfg.dt_myr as f64, cfl_cfg);
        cfl.update(&res);
    }
    // Transforms raw deltas
    let mut tr_tmp = vec![0.0f32; w.grid.cells];
    let _ = engine::transforms::apply_transforms(
        &w.grid,
        &w.boundaries,
        &w.plates.plate_id,
        &w.v_en,
        &mut tr_tmp,
        engine::transforms::TransformParams {
            tau_open_m_per_yr: 0.005,
            min_tangential_m_per_yr: 0.015,
            max_normal_m_per_yr: 0.001,
            basin_half_width_km: 10.0,
            ridge_like_uplift_m: -3.0,
            basin_deepen_m: 6.0,
        },
        cfg.dt_myr as f64,
    );
    let width_tr_m = (50.0f64).clamp(20.0, 80.0) * 1000.0;
    for &raw in &tr_tmp {
        let res = engine::cfl::limit(raw as f64, width_tr_m, cfg.dt_myr as f64, cfl_cfg);
        cfl.update(&res);
    }

    // Basic sanity and determinism: values finite and repeatable across a re-run
    assert!(land_frac.is_finite());
    assert!(n_r + n_s + n_t > 0);
    if let Some((_cap_comp, _cap_thc, _rs, _ru)) = caps { /* present and parsed */ }
    assert!(cfl.max_cfl.is_finite() && cfl.mean_cfl.is_finite());

    // Re-run and compare within tight tolerance
    let mut w2 = World::new(f, plates, seed);
    let mut elev2 = vec![0.0f32; w2.grid.cells];
    let mut eta2: f32 = 0.0;
    let surf2 = pipeline::SurfaceFields { elevation_m: &mut elev2, eta_m: &mut eta2 };
    pipeline::step_full(&mut w2, surf2, cfg);
    // boundary counts
    let mut n_r2 = 0u32;
    let mut n_s2 = 0u32;
    let mut n_t2 = 0u32;
    for &(_u, _v, cls) in &w2.boundaries.edges {
        match cls {
            1 => n_r2 += 1,
            2 => n_s2 += 1,
            3 => n_t2 += 1,
            _ => {}
        }
    }
    assert_eq!(n_r, n_r2);
    assert_eq!(n_s, n_s2);
    assert_eq!(n_t, n_t2);
    // land fraction
    let mut land2 = 0.0f64;
    let mut tot2 = 0.0f64;
    for (a, &z) in w2.area_m2.iter().zip(elev2.iter()) {
        tot2 += *a as f64;
        if z > 0.0 {
            land2 += *a as f64;
        }
    }
    let land_frac2 = if tot2 > 0.0 { land2 / tot2 } else { 0.0 };
    assert!((land_frac - land_frac2).abs() <= 1e-9);
}
