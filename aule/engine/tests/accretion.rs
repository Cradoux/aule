use engine as crate_engine;

#[test]
fn oc_no_change_when_no_convergence_or_both_oceanic() {
    let f: u32 = 8;
    let mut world = crate_engine::world::World::new(f, 2, 123);
    // Make both sides oceanic: C near zero
    for ci in &mut world.c {
        *ci = 0.0;
    }
    // Subduction/masks
    let sub = crate_engine::subduction::apply_subduction(
        &world.grid,
        &world.boundaries,
        &world.plates.plate_id,
        &world.age_myr,
        &world.v_en,
        &mut world.depth_m,
        crate_engine::subduction::SubductionParams {
            tau_conv_m_per_yr: 0.005,
            trench_half_width_km: 50.0,
            arc_offset_km: 150.0,
            arc_half_width_km: 30.0,
            backarc_width_km: 150.0,
            trench_deepen_m: 0.0,
            arc_uplift_m: 0.0,
            backarc_uplift_m: 0.0,
            rollback_offset_m: 0.0,
            rollback_rate_km_per_myr: 0.0,
            backarc_extension_mode: false,
            backarc_extension_deepen_m: 0.0,
        },
    );
    let c0 = world.c.clone();
    let th0 = world.th_c_m.clone();
    let p = crate_engine::accretion::AccretionParams {
        k_arc: 0.05,
        gamma_obliquity: 1.0,
        beta_arc: 1.0,
        alpha_arc: 0.02,
        alpha_forearc: 0.01,
        c_min_continent: 0.6,
        thc_min_m: 20_000.0,
        thc_max_m: 70_000.0,
        enable_docking: false,
        c_terrane_min: 0.5,
        d_dock_km: 150.0,
        vn_min_m_per_yr: 0.005,
        tau_dock: 0.02,
        couple_flexure: false,
    };
    let s = crate_engine::accretion::apply_oc_accretion(
        &world.grid,
        &sub.masks,
        &world.boundaries,
        &crate_engine::plates::velocity_field_m_per_yr(
            &world.grid,
            &world.plates,
            &world.plates.plate_id,
        ),
        &mut world.c,
        &mut world.th_c_m,
        &mut world.depth_m,
        &world.area_m2,
        &p,
        1.0,
    );
    assert_eq!(s.edges_oc, 0);
    assert_eq!(s.cells_arc + s.cells_forearc, 0);
    assert!(world.c.iter().zip(c0.iter()).all(|(a, b)| (*a - *b).abs() < 1e-6));
    assert!(world.th_c_m.iter().zip(th0.iter()).all(|(a, b)| (*a - *b).abs() < 1e-6));
}
