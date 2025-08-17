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
fn offset_solves_ocean_fraction_small() {
    // Mixed depths and areas, check ocean fraction within 1e-4
    let depth =
        vec![200.0f32, -50.0, 0.0, 300.0, -10.0, 800.0, 20.0, -5.0, 0.0, 100.0, -200.0, 50.0];
    let area: Vec<f32> = (0..depth.len()).map(|i| 1.0 + (i as f32) * 0.1).collect();
    let total: f64 = area.iter().map(|&a| a as f64).sum();
    let target_ocean = 0.65f32; // 65% ocean
    let off = isostasy::solve_offset_for_ocean_area_fraction(&depth, &area, target_ocean, 1e-4, 64);
    let mut ocean = 0.0f64;
    for (d, a) in depth.iter().zip(area.iter()) {
        if (*d as f64 + off) > 0.0 {
            ocean += *a as f64;
        }
    }
    let frac = if total > 0.0 { ocean / total } else { 0.0 };
    assert!((frac - target_ocean as f64).abs() <= 1e-4);
}

#[test]
fn amplitude_solves_land_fraction_and_is_deterministic() {
    // Synthesize a simple continent template with two caps (handmade)
    let n = 100usize;
    let mut tpl = vec![0.0f32; n];
    for t in tpl.iter_mut().take(40).skip(20) {
        *t = 1.0;
    }
    for t in tpl.iter_mut().take(75).skip(60) {
        *t = 0.7;
    }
    let base: Vec<f32> = (0..n).map(|i| if i % 3 == 0 { 3000.0 } else { 2000.0 }).collect();
    let area: Vec<f32> = vec![2.0; n];
    let target_land = 0.30f32;
    let (amp_a, off_a) = isostasy::solve_amplitude_for_land_fraction(
        &tpl,
        &base,
        &area,
        target_land,
        0.0,
        6000.0,
        2e-2,
        1e-4,
        48,
    );
    let (amp_b, off_b) = isostasy::solve_amplitude_for_land_fraction(
        &tpl,
        &base,
        &area,
        target_land,
        0.0,
        6000.0,
        2e-2,
        1e-4,
        48,
    );
    assert_eq!(amp_a, amp_b);
    assert_eq!(off_a, off_b);
    // Apply and check achieved land fraction
    let mut tmp = base.clone();
    for (d, t) in tmp.iter_mut().zip(tpl.iter()) {
        *d = *d - (amp_a as f32) * *t + (off_a as f32);
    }
    let total: f64 = area.iter().map(|&a| a as f64).sum();
    let mut land = 0.0f64;
    for (d, a) in tmp.iter().zip(area.iter()) {
        if *d <= 0.0 {
            land += *a as f64;
        }
    }
    let land_frac = if total > 0.0 { land / total } else { 0.0 };
    assert!((land_frac - target_land as f64).abs() <= 0.02);
}

#[test]
fn solver_hits_target_within_tol() {
    // Synthetic: depths 0..999, area=2.0 each
    let n = 1000usize;
    let mut depth = vec![0.0f32; n];
    let area = vec![2.0f32; n];
    for (i, d) in depth.iter_mut().enumerate().take(n) {
        *d = (i as f32) * 1.0;
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
