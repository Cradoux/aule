use engine::{age, boundaries::Boundaries, grid::Grid, plates::Plates, subduction};

fn mean(values: &[f32]) -> f32 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f32>() / values.len() as f32
    }
}

#[test]
fn determinism_and_idempotence() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 2, 123);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    let aout = age::compute_age_and_bathymetry(
        &g,
        &b,
        &plates.plate_id,
        &plates.vel_en,
        age::AgeParams::default(),
    );
    let mut depth = aout.depth_m.clone();
    let p = subduction::SubductionParams {
        tau_conv_m_per_yr: 0.005,
        trench_half_width_km: 30.0,
        arc_offset_km: 150.0,
        arc_half_width_km: 30.0,
        backarc_width_km: 150.0,
        trench_deepen_m: 100.0,
        arc_uplift_m: -20.0,
        backarc_uplift_m: -10.0,
    };
    let r1 = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth,
        p,
    );
    let depth_after_once = depth.clone();
    let r2 = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth,
        p,
    );
    assert_eq!(r1.stats.trench_cells, r2.stats.trench_cells);
    assert_eq!(r1.stats.arc_cells, r2.stats.arc_cells);
    assert_eq!(r1.stats.backarc_cells, r2.stats.backarc_cells);
    assert_eq!(depth_after_once, depth);
}

#[test]
fn bathy_signs() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 2, 321);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    let aout = age::compute_age_and_bathymetry(
        &g,
        &b,
        &plates.plate_id,
        &plates.vel_en,
        age::AgeParams::default(),
    );
    let mut depth = aout.depth_m.clone();
    let p = subduction::SubductionParams {
        tau_conv_m_per_yr: 0.005,
        trench_half_width_km: 30.0,
        arc_offset_km: 150.0,
        arc_half_width_km: 30.0,
        backarc_width_km: 150.0,
        trench_deepen_m: 500.0,
        arc_uplift_m: -100.0,
        backarc_uplift_m: -50.0,
    };
    let base = depth.clone();
    let r = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth,
        p,
    );
    if r.stats.trench_cells > 0 {
        let mut trench_vals: Vec<f32> = Vec::new();
        for (i, &m) in r.masks.trench.iter().enumerate() {
            if m {
                trench_vals.push(depth[i] - base[i]);
            }
        }
        assert!(mean(&trench_vals) > 0.0);
    }
    if r.stats.arc_cells > 0 {
        let mut arc_vals: Vec<f32> = Vec::new();
        for (i, &m) in r.masks.arc.iter().enumerate() {
            if m {
                arc_vals.push(depth[i] - base[i]);
            }
        }
        assert!(mean(&arc_vals) < 0.0);
    }
    if r.stats.backarc_cells > 0 {
        let mut bavals: Vec<f32> = Vec::new();
        for (i, &m) in r.masks.backarc.iter().enumerate() {
            if m {
                bavals.push(depth[i] - base[i]);
            }
        }
        assert!(mean(&bavals) < 0.0);
    }
}
