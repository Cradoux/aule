use engine as crate_engine;

#[test]
fn cc_gates_by_cmin() {
    let f: u32 = 8;
    let mut world = crate_engine::world::World::new(f, 2, 42);
    // Make a fake convergent boundary by forcing small opening threshold in classifier is already present
    // Ensure one side oceanic
    for i in 0..world.grid.cells {
        world.c[i] = if i % 2 == 0 { 0.7 } else { 0.2 };
    }
    // Build velocities from plates
    let vel3 = crate_engine::plates::velocity_field_m_per_yr(
        &world.grid,
        &world.plates,
        &world.plates.plate_id,
    );
    // Save baseline
    let th_before = world.th_c_m.clone();
    let p = crate_engine::orogeny::OrogenyParams {
        c_min: 0.6,
        w_core_km: 100.0,
        w_taper_km: 100.0,
        k_thick: 0.1,
        beta_uplift: 1.0,
        gamma_obliquity: 1.0,
        couple_flexure: false,
    };
    let _ = crate_engine::orogeny::apply_cc_orogeny(
        &world.grid,
        &world.boundaries,
        &world.plates.plate_id,
        &vel3,
        &world.c,
        &world.area_m2,
        &mut world.th_c_m,
        &mut world.depth_m,
        &p,
        1.0,
    );
    // Oceanic side suppresses effect => most cells unchanged in this crude test
    let mut changed = 0usize;
    for (i, &th0) in th_before.iter().enumerate().take(world.grid.cells) {
        if (world.th_c_m[i] - th0).abs() > 1e-3 {
            changed += 1;
        }
    }
    // still allow a few due to coarse constraints
    assert!(changed < world.grid.cells / 2);
}
