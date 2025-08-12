use engine::{continent, grid::Grid};

#[test]
fn determinism_template() {
    let g = Grid::new(16);
    let p = continent::ContinentParams {
        seed: 42,
        n_continents: 3,
        mean_radius_km: 2200.0,
        falloff_km: 600.0,
        plateau_uplift_m: -2500.0,
        target_land_fraction: None,
    };
    let c0 = continent::build_continents(&g, p);
    let c1 = continent::build_continents(&g, p);
    assert_eq!(c0.uplift_template_m, c1.uplift_template_m);
}

#[test]
fn amplitude_targets_land_fraction() {
    let g = Grid::new(16);
    // synthetic depths: radial with some spread
    let mut depth = vec![3000.0f32; g.cells];
    // areas in m^2 (unit-sphere areas scaled by Earth sphere)
    const R: f64 = 6_371_000.0;
    let scale = 4.0 * std::f64::consts::PI * R * R;
    let area_m2: Vec<f32> = g.area.iter().map(|&a| (a as f64 * scale) as f32).collect();
    let p = continent::ContinentParams {
        seed: 7,
        n_continents: 3,
        mean_radius_km: 2200.0,
        falloff_km: 600.0,
        plateau_uplift_m: -2500.0,
        target_land_fraction: None,
    };
    let c = continent::build_continents(&g, p);
    // With these cap sizes the reachable union is ~8–10%
    let target = 0.08;
    let amp = continent::solve_amplitude_for_target_land_fraction(
        &depth,
        &c.uplift_template_m,
        &area_m2,
        target,
        1e-3,
        128,
    );
    let (_mask, frac) =
        continent::apply_continents(&mut depth, &c.uplift_template_m, amp, &area_m2);
    assert!((frac - target).abs() <= 0.02, "frac={} target={}", frac, target);
}
