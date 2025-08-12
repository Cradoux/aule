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
