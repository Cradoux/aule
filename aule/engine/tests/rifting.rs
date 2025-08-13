use engine as crate_engine;

#[test]
fn no_divergence_no_rifting() {
    let f: u32 = 8;
    let mut world = crate_engine::world::World::new(f, 2, 7);
    // Make everything oceanic to gate out rifting
    for ci in &mut world.c {
        *ci = 0.0;
    }
    let c0 = world.c.clone();
    let th0 = world.th_c_m.clone();
    let d0 = world.depth_m.clone();
    let stats = crate_engine::rifting::apply_rifting(
        &world.grid,
        &world.boundaries,
        &world.plates.plate_id,
        &crate_engine::plates::velocity_field_m_per_yr(
            &world.grid,
            &world.plates,
            &world.plates.plate_id,
        ),
        &mut world.c,
        &mut world.th_c_m,
        &mut world.age_myr,
        &mut world.depth_m,
        &world.area_m2,
        &crate_engine::rifting::RiftingParams {
            c_rift_min: 0.6,
            v_open_min_m_per_yr: 0.001,
            w_core_km: 60.0,
            w_taper_km: 250.0,
            k_thin: 0.15,
            alpha_subs: 1.0,
            ocean_thresh: 0.15,
            k_c_oceanize: 0.05,
            reset_age_on_core: true,
            enable_shoulder: false,
            w_bulge_km: 120.0,
            beta_shoulder: 0.3,
            couple_flexure: false,
            thc_min_m: 20_000.0,
            thc_max_m: 70_000.0,
        },
        1.0,
    );
    assert_eq!(stats.edges_rift, 0);
    assert_eq!(stats.cells_core + stats.cells_margin, 0);
    assert!(world.c.iter().zip(c0.iter()).all(|(a, b)| (*a - *b).abs() < 1e-6));
    assert!(world.th_c_m.iter().zip(th0.iter()).all(|(a, b)| (*a - *b).abs() < 1e-6));
    assert!(world.depth_m.iter().zip(d0.iter()).all(|(a, b)| (*a - *b).abs() < 1e-6));
}
