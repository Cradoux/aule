use engine::{
    flexure_loads::{assemble_load_from_depth, LoadParams},
    grid::Grid,
};

#[test]
fn units_sanity_and_continuity() {
    let g = Grid::new(8);
    let p = LoadParams { rho_w: 1030.0, rho_c: 2900.0, g: 9.81, sea_level_m: 0.0 };
    let mut depth = vec![0.0f32; g.cells];
    // 3 km ocean everywhere
    for d in &mut depth {
        *d = 3000.0;
    }
    let f = assemble_load_from_depth(&g, &depth, &p);
    // ρ_w g h = 1030*9.81*3000 ≈ 3.03e7 N/m^2
    let expected = 1030.0 * 9.81 * 3000.0;
    let mean: f32 = f.iter().copied().sum::<f32>() / f.len().max(1) as f32;
    let rel = (mean - expected).abs() / expected;
    println!("mean={} expected={} rel={}", mean, expected, rel);
    assert!(rel < 1e-4);

    // continuity around 0: left/right limits match coefficients at ±ε
    let eps = 1e-3f32;
    let f_neg = assemble_load_from_depth(&g, &vec![-eps; g.cells], &p)[0];
    let f_pos = assemble_load_from_depth(&g, &vec![eps; g.cells], &p)[0];
    let slope_left = (f_neg - 0.0) / eps; // ≈ rho_c g
    let slope_right = (f_pos - 0.0) / eps; // ≈ rho_w g
    assert!((slope_left - p.rho_c * p.g).abs() / (p.rho_c * p.g) < 1e-6);
    assert!((slope_right - p.rho_w * p.g).abs() / (p.rho_w * p.g) < 1e-6);
}

#[test]
fn determinism() {
    let g = Grid::new(4);
    let mut depth = vec![0.0f32; g.cells];
    for (i, d) in depth.iter_mut().enumerate() {
        *d = if (i % 2) == 0 { 1000.0 } else { -500.0 };
    }
    let p = LoadParams { rho_w: 1000.0, rho_c: 2800.0, g: 9.81, sea_level_m: 0.0 };
    let a = assemble_load_from_depth(&g, &depth, &p);
    let b = assemble_load_from_depth(&g, &depth, &p);
    assert_eq!(a, b);
}
