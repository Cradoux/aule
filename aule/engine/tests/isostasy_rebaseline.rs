use engine as crate_engine;

#[test]
fn compute_ref_deterministic() {
    // Fixed inputs
    let depth = vec![1000.0f32, -200.0, 300.0, 0.0, 50.0];
    let area = vec![1.0e6f32; depth.len()];
    let eta = 0.0f32;
    let r1 = crate_engine::isostasy::compute_ref(&depth, &area, eta);
    let r2 = crate_engine::isostasy::compute_ref(&depth, &area, eta);
    assert_eq!(r1.ocean_area_m2, r2.ocean_area_m2);
    assert_eq!(r1.volume_m3, r2.volume_m3);
}

#[test]
fn zero_offset_after_rebaseline() {
    // Construct a simple world
    let f: u32 = 8;
    let mut world = crate_engine::world::World::new(f, 2, 123);
    // Synthetic depths & areas
    for d in &mut world.depth_m {
        *d = 1000.0;
    }
    // Rebaseline and then solve for offset; expect near zero
    let area = world.area_m2.clone();
    let _ = crate_engine::isostasy::rebaseline(&mut world, &area);
    let r = world.sea_level_ref.unwrap();
    let off = crate_engine::isostasy::solve_offset_for_volume(
        &world.depth_m,
        &world.area_m2,
        r.volume_m3,
        1e6,
        64,
    );
    assert!(off.abs() < 1e-3, "offset not ~0 after rebaseline: {}", off);
}

#[test]
fn land_fraction_preserved_one_solve_after_rebaseline() {
    // Build world with mix of land and ocean
    let f: u32 = 8;
    let mut world = crate_engine::world::World::new(f, 2, 123);
    for (i, d) in world.depth_m.iter_mut().enumerate() {
        *d = if i % 3 == 0 { -100.0 } else { 500.0 };
    }
    let land_frac_before = {
        let mut land_area = 0.0f64;
        let mut tot = 0.0f64;
        for i in 0..world.depth_m.len() {
            tot += world.area_m2[i] as f64;
            if world.depth_m[i] <= 0.0 {
                land_area += world.area_m2[i] as f64;
            }
        }
        if tot > 0.0 {
            land_area / tot
        } else {
            0.0
        }
    };
    let area = world.area_m2.clone();
    let _ = crate_engine::isostasy::rebaseline(&mut world, &area);
    let r = world.sea_level_ref.unwrap();
    let off = crate_engine::isostasy::solve_offset_for_volume(
        &world.depth_m,
        &world.area_m2,
        r.volume_m3,
        1e6,
        64,
    );
    crate_engine::isostasy::apply_sea_level_offset(&mut world.depth_m, off);
    let land_frac_after = {
        let mut land_area = 0.0f64;
        let mut tot = 0.0f64;
        for i in 0..world.depth_m.len() {
            tot += world.area_m2[i] as f64;
            if world.depth_m[i] <= 0.0 {
                land_area += world.area_m2[i] as f64;
            }
        }
        if tot > 0.0 {
            land_area / tot
        } else {
            0.0
        }
    };
    let diff = (land_frac_before - land_frac_after).abs();
    assert!(
        diff < 1e-6,
        "land fraction changed: {} -> {} (diff={})",
        land_frac_before,
        land_frac_after,
        diff
    );
}
