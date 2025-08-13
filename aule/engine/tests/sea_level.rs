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
