use engine::world::{SimplePreset, World};

#[test]
fn simple_generate_variance_and_determinism() {
    let f = 16u32;
    let seed = 1_234_567u64;
    let mut w = World::new(f, 8, seed);
    let preset = SimplePreset { plates: 8, continents_n: 3, target_land_frac: 0.30 };
    let _rep = w.generate_simple(&preset, seed);
    // Variance check (rough): ensure depth is non-flat
    let mut mean = 0.0f64;
    let n = w.depth_m.len() as f64;
    for &d in &w.depth_m {
        mean += d as f64;
    }
    mean /= n.max(1.0);
    let mut var = 0.0f64;
    for &d in &w.depth_m {
        let x = (d as f64) - mean;
        var += x * x;
    }
    assert!(var > 1e5, "variance too small: {}", var);

    // Determinism: second run with same seed & F yields identical depth
    let mut w2 = World::new(f, 8, seed);
    let _rep2 = w2.generate_simple(&preset, seed);
    assert_eq!(w.depth_m.len(), w2.depth_m.len());
    for i in 0..w.depth_m.len() {
        assert!((w.depth_m[i] - w2.depth_m[i]).abs() < 1e-6);
    }
}
