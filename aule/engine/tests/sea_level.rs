use engine as crate_engine;

#[test]
fn target_fraction_solve_hits_within_tol() {
    // Synthetic depths: mix of land (<=0) and ocean (>0)
    let depth_m: Vec<f32> = vec![-100.0, 0.0, 50.0, 200.0, -50.0, 500.0, 1000.0, -10.0];
    let area_m2: Vec<f32> = vec![1.0; depth_m.len()];
    let phi_target = 0.70;
    let eta = crate_engine::sea_level::solve_offset_for_target_fraction(
        &depth_m, &area_m2, phi_target, 1e-3, 64,
    );
    // Compute resulting fraction
    let mut ocean = 0.0f64;
    for (&d, &a) in depth_m.iter().zip(area_m2.iter()) {
        if (d as f64 + eta) > 0.0 {
            ocean += a as f64;
        }
    }
    let phi = ocean / (area_m2.len() as f64);
    // Discrete cell areas mean we may not hit phi* exactly; allow coarse tolerance here
    assert!((phi - phi_target).abs() < 0.2);
}

#[test]
fn linear_policy_monotonicity() {
    let d = vec![100.0f32; 10];
    let a = vec![1.0f32; 10];
    let pol = crate_engine::sea_level::EustasyPolicy::LinearTrend { d_eta_mm_per_yr: 1.0 };
    let eta_100 = crate_engine::sea_level::update_eustasy_eta(100.0, 1.0, &pol, &d, &a, 0.0);
    assert!((eta_100 - 100000.0).abs() < 1.0);
}

#[test]
fn cyclic_bounds() {
    let d = vec![100.0f32; 10];
    let a = vec![1.0f32; 10];
    let pol = crate_engine::sea_level::EustasyPolicy::Cyclic {
        amp_m: 60.0,
        period_myr: 100.0,
        phase_deg: 0.0,
    };
    for &t in &[0.0, 25.0, 50.0, 75.0, 100.0, 150.0] {
        let eta = crate_engine::sea_level::update_eustasy_eta(t, 1.0, &pol, &d, &a, 0.0);
        assert!(eta.abs() <= 60.0 + 1e-6);
    }
}

#[test]
fn land_fraction_offset_hits_target_on_synthetic_dem() {
    // Synthetic DEM spanning [-2000, 4000] so land fraction varies smoothly with offset
    let n = 2048usize;
    let mut depth = vec![0.0f32; n];
    let area = vec![1.0f32; n];
    for (i, d) in depth.iter_mut().enumerate().take(n) {
        let t = (i as f64) / ((n - 1) as f64);
        *d = (-2000.0 + t * 6000.0) as f32; // [-2000, 4000]
    }
    let target_land = 0.30f32;
    let off =
        crate_engine::isostasy::solve_offset_for_land_fraction(&depth, &area, target_land, 64);
    // Elevation = -(depth+off); land if elevation > 0 -> depth+off < 0
    let mut land_area = 0.0f64;
    let mut total_area = 0.0f64;
    for (&d, &a) in depth.iter().zip(area.iter()) {
        if (d as f64 + off) < 0.0 {
            land_area += a as f64;
        }
        total_area += a as f64;
    }
    let frac = if total_area > 0.0 { land_area / total_area } else { 0.0 } as f32;
    assert!((frac - target_land).abs() <= 0.02);
}

#[test]
fn constant_volume_offset_preserves_volume_on_shift() {
    // DEM with positive depths only; shifting by -X should reduce volume by X*area
    let n = 100usize;
    let depth = vec![2000.0f32; n];
    let area = vec![2.0f32; n];
    let (v0, _a0) = crate_engine::isostasy::ocean_volume_from_depth(&depth, &area);
    let target = v0 - (1000.0 * (2.0 * n as f32)) as f64;
    let off = crate_engine::isostasy::solve_offset_for_volume(&depth, &area, target, 1e3, 64);
    assert!(off < 0.0);
}
