use engine::{age, boundaries::Boundaries, grid::Grid, plates::Plates};

#[test]
fn determinism_and_nonnegativity() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 4, 1234);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    let p = age::AgeParams { v_floor_m_per_yr: 0.005 };
    let out1 = age::compute_age_and_bathymetry(&g, &b, &plates.plate_id, &plates.vel_en, p);
    let out2 = age::compute_age_and_bathymetry(&g, &b, &plates.plate_id, &plates.vel_en, p);
    assert_eq!(out1.age_myr.len(), g.cells);
    assert_eq!(out1.depth_m.len(), g.cells);
    for i in 0..g.cells {
        assert!(out1.age_myr[i] >= 0.0);
        assert!(out1.depth_m[i] >= 0.0);
    }
    assert_eq!(out1.age_myr, out2.age_myr);
    assert_eq!(out1.depth_m, out2.depth_m);
}

#[test]
fn depth_from_age_monotonic_and_anchors() {
    let d0 = 2600.0;
    let a = 350.0;
    let b = 0.0;
    // Anchors
    let d_zero = age::depth_from_age(0.0, d0, a, b);
    assert!((d_zero - d0).abs() < 1e-6);
    // Monotonic increasing with age
    let d_10 = age::depth_from_age(10.0, d0, a, b);
    let d_20 = age::depth_from_age(20.0, d0, a, b);
    assert!(d_20 > d_10);
    // Rough anchor around 100 Myr
    let d_100 = age::depth_from_age(100.0, d0, a, b);
    let expected = d0 + a * (100.0_f64).sqrt() + b * 100.0;
    assert!((d_100 - expected).abs() < 1e-6);
}
