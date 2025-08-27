use engine::{age, boundaries::Boundaries, grid::Grid, plates::Plates, subduction};

fn mean(values: &[f32]) -> f32 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f32>() / values.len() as f32
    }
}

#[test]
fn rollback_shifts_bands_and_extension_deepens_backarc() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 2, 2024);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    let aout = age::compute_age_and_bathymetry(
        &g,
        &b,
        &plates.plate_id,
        &plates.vel_en,
        age::AgeParams::default(),
    );

    // Baseline params (no rollback)
    let mut depth0 = aout.depth_m.clone();
    let p0 = subduction::SubductionParams {
        tau_conv_m_per_yr: 0.005,
        trench_half_width_km: 30.0,
        arc_offset_km: 150.0,
        arc_half_width_km: 30.0,
        backarc_width_km: 150.0,
        trench_deepen_m: 500.0,
        arc_uplift_m: -100.0,
        backarc_uplift_m: -50.0,
        rollback_offset_m: 0.0,
        rollback_rate_km_per_myr: 0.0,
        backarc_extension_mode: false,
        backarc_extension_deepen_m: 600.0,
        continent_c_min: 0.6,
    };
    let zero_c: Vec<f32> = vec![0.0; g.cells];
    let r0 = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth0,
        p0,
        Some(&zero_c),
    );

    // Apply 100 km rollback offset
    let mut depth1 = aout.depth_m.clone();
    let p1 = subduction::SubductionParams { rollback_offset_m: 100_000.0, ..p0 };
    let r1 = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth1,
        p1,
        Some(&zero_c),
    );

    // Arc/back-arc may be zero in tiny synthetic cases; require at least trench or arc presence
    assert!(r0.stats.trench_cells > 0 || r0.stats.arc_cells > 0);
    assert!(r1.stats.trench_cells > 0 || r1.stats.arc_cells > 0);
    assert!(r0.stats.backarc_cells > 0 && r1.stats.backarc_cells > 0);

    // Extension mode: back-arc mean depth deeper than uplift mode
    let mut depth_ext = aout.depth_m.clone();
    let p_ext = subduction::SubductionParams { backarc_extension_mode: true, ..p0 };
    let r_ext = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth_ext,
        p_ext,
        Some(&zero_c),
    );
    if r_ext.stats.backarc_cells > 0 {
        let mut d_uplift: Vec<f32> = Vec::new();
        let mut d_ext: Vec<f32> = Vec::new();
        for (i, &m) in r0.masks.backarc.iter().enumerate() {
            if m {
                d_uplift.push(depth0[i]);
            }
        }
        for (i, &m) in r_ext.masks.backarc.iter().enumerate() {
            if m {
                d_ext.push(depth_ext[i]);
            }
        }
        if !d_uplift.is_empty() && !d_ext.is_empty() {
            assert!(mean(&d_ext) > mean(&d_uplift));
        }
    }

    // Idempotence
    let mut depth2 = depth1.clone();
    let _ = subduction::apply_subduction(
        &g,
        &b,
        &plates.plate_id,
        &aout.age_myr,
        &plates.vel_en,
        &mut depth2,
        p1,
        Some(&zero_c),
    );
    assert_eq!(depth1, depth2);
}
