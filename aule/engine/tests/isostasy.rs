use engine::isostasy;

#[test]
fn volume_monotonic_in_offset() {
    let depth = vec![0.0f32, 1000.0, 2000.0, 0.0, 3000.0];
    let area = vec![1.0f32; depth.len()];
    let v1 = isostasy::solve_offset_for_volume(&depth, &area, 2_000.0, 1e-9, 128);
    let v2 = isostasy::solve_offset_for_volume(&depth, &area, 3_000.0, 1e-9, 128);
    assert!(v2 > v1);
}

#[test]
fn solver_hits_target_within_tol() {
    // Synthetic: depths 0..999, area=2.0 each
    let n = 1000usize;
    let mut depth = vec![0.0f32; n];
    let area = vec![2.0f32; n];
    for i in 0..n {
        depth[i] = (i as f32) * 1.0;
    }
    let (v0, _a0) = isostasy::ocean_volume_from_depth(&depth, &area);
    // Target: add 100 m average depth → additional volume = 100 * sum(area)
    let sum_area: f64 = area.iter().map(|a| *a as f64).sum();
    let target = v0 + 100.0 * sum_area;
    let off = isostasy::solve_offset_for_volume(&depth, &area, target, 1e-6 * target, 128);
    let after = volume_with_offset_ref(&depth, &area, off);
    let err = (after - target).abs();
    assert!(err <= 1e-6 * target, "err={} target={}", err, target);
}

#[test]
fn idempotence_apply_offset() {
    let mut depth = vec![100.0f32, 0.0, 50.0];
    isostasy::apply_sea_level_offset(&mut depth, 20.0);
    let once = depth.clone();
    isostasy::apply_sea_level_offset(&mut depth, 0.0);
    assert_eq!(depth, once);
}

#[test]
fn determinism_same_inputs_same_offset() {
    let depth = vec![0.0f32, 1000.0, 2000.0, 100.0];
    let area = vec![1.0f32, 2.0, 3.0, 4.0];
    let target = 10_000.0f64;
    let a = isostasy::solve_offset_for_volume(&depth, &area, target, 1e-9, 128);
    let b = isostasy::solve_offset_for_volume(&depth, &area, target, 1e-9, 128);
    assert_eq!(a, b);
}

fn volume_with_offset_ref(depth: &[f32], area: &[f32], off: f64) -> f64 {
    let mut v = 0.0f64;
    for i in 0..depth.len() {
        let d = (depth[i] as f64) + off;
        if d > 0.0 {
            v += d * (area[i] as f64);
        }
    }
    v
}
