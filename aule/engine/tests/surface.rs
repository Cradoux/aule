//! Tests for surface processes (T-450).

use engine::surface;
use engine::world::World;

#[test]
fn zero_toggles_no_change() {
    let mut world = World::new(8, 4, 1234);
    // Flat land small bump
    for d in &mut world.depth_m {
        *d = 0.0;
    }
    world.depth_m[0] = -100.0;
    let before = world.depth_m.clone();
    let stats = surface::apply_surface_processes(
        &world.grid,
        &world.c,
        &mut world.depth_m,
        &mut world.sediment_m,
        &world.area_m2,
        &surface::SurfaceParams { k_stream: 0.0, k_diff: 0.0, k_tr: 0.0, ..Default::default() },
        1.0,
    );
    assert!(stats.eroded_m3.abs() < 1e-12);
    assert!(stats.deposited_m3.abs() < 1e-12);
    assert_eq!(world.depth_m, before);
}
