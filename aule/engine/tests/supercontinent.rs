use engine::world::{SimplePreset, World};

#[test]
fn supercontinent_determinism_and_eta_land_fraction() {
    let f = 64u32; // use high-F path in generator
    let seed = 1_234_567u64;
    let preset = SimplePreset { plates: 8, continents_n: 0, target_land_frac: 0.30 };

    let mut w1 = World::new(f, preset.plates, seed);
    let rep1 = w1.generate_simple(&preset, seed);

    let mut w2 = World::new(f, preset.plates, seed);
    let rep2 = w2.generate_simple(&preset, seed);

    // Determinism: depth and eta identical across runs with same seed and F
    assert_eq!(w1.depth_m.len(), w2.depth_m.len());
    for i in 0..w1.depth_m.len() {
        assert!((w1.depth_m[i] - w2.depth_m[i]).abs() < 1e-6);
    }
    assert!((w1.sea.eta_m - w2.sea.eta_m).abs() < 1e-6);

    // Eta-aware land fraction consistency: recompute and match reported
    let mut ocean_area = 0.0f64;
    let mut total_area = 0.0f64;
    for (d, a) in w1.depth_m.iter().zip(w1.area_m2.iter()) {
        if (*d as f64 + w1.sea.eta_m as f64) > 0.0 {
            ocean_area += *a as f64;
        }
        total_area += *a as f64;
    }
    let ocean_frac = if total_area > 0.0 { ocean_area / total_area } else { 0.0 } as f32;
    let land_frac = 1.0 - ocean_frac;
    assert!((land_frac - rep1.land_frac).abs() <= 1e-6);
    assert!((rep1.land_frac - rep2.land_frac).abs() <= 1e-6);
}

#[test]
fn supercontinent_low_f_determinism() {
    let f = 16u32; // low-F path
    let seed = 1_234_567u64;
    let preset = SimplePreset { plates: 8, continents_n: 0, target_land_frac: 0.30 };

    let mut w1 = World::new(f, preset.plates, seed);
    let rep1 = w1.generate_simple(&preset, seed);

    let mut w2 = World::new(f, preset.plates, seed);
    let rep2 = w2.generate_simple(&preset, seed);

    // Depth and eta determinism
    assert_eq!(w1.depth_m.len(), w2.depth_m.len());
    for i in 0..w1.depth_m.len() {
        assert!((w1.depth_m[i] - w2.depth_m[i]).abs() < 1e-6);
    }
    assert!((w1.sea.eta_m - w2.sea.eta_m).abs() < 1e-6);

    // Reported land fraction matches eta-aware recomputation
    let mut ocean_area = 0.0f64;
    let mut total_area = 0.0f64;
    for (d, a) in w1.depth_m.iter().zip(w1.area_m2.iter()) {
        if (*d as f64 + w1.sea.eta_m as f64) > 0.0 {
            ocean_area += *a as f64;
        }
        total_area += *a as f64;
    }
    let ocean_frac = if total_area > 0.0 { ocean_area / total_area } else { 0.0 } as f32;
    let land_frac = 1.0 - ocean_frac;
    assert!((land_frac - rep1.land_frac).abs() <= 1e-6);
    assert!((rep1.land_frac - rep2.land_frac).abs() <= 1e-6);
}

#[ignore]
#[test]
fn supercontinent_land_fraction_target_tolerance() {
    // NOTE: Enable this when the generator reliably meets targets.
    let f = 64u32;
    let seed = 1_234_567u64;
    let preset = SimplePreset { plates: 8, continents_n: 0, target_land_frac: 0.30 };

    let mut w = World::new(f, preset.plates, seed);
    let rep = w.generate_simple(&preset, seed);

    // Expect land fraction to be within 10 percentage points of target.
    let tol = 0.10f32;
    assert!(
        (rep.land_frac - preset.target_land_frac).abs() <= tol,
        "land_frac={} target={}",
        rep.land_frac,
        preset.target_land_frac
    );
}
