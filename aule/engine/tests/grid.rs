use engine::grid::{cache::GridError, Grid};
use std::{
    env, fs,
    time::{SystemTime, UNIX_EPOCH},
};

#[test]
fn grid_constructs() {
    let g = Grid::new(0);
    assert_eq!(g.cells, 12);
    assert_eq!(g.pos_xyz.len(), g.cells);
    assert_eq!(g.latlon.len(), g.cells);
    assert_eq!(g.area.len(), g.cells);
    assert_eq!(g.n1.len(), g.cells);
    assert_eq!(g.n2.len(), g.cells);
}

fn rel_err(a: f64, b: f64) -> f64 {
    (a - b).abs() / b
}

#[test]
fn area_partition_f16_f32_and_counts() {
    for &f in &[16u32, 32u32] {
        let g = Grid::new(f);
        // Total area ~ 4π
        let sum_area: f64 = g.area.iter().map(|a| *a as f64).sum();
        let sphere = 4.0 * std::f64::consts::PI;
        assert!(rel_err(sum_area, sphere) < 1e-6, "total area rel err too large for F={f}");

        // Degree counts: exactly 12 pentagons (5), rest hexagons (6)
        let pent = g.n1.iter().filter(|n| n.len() == 5).count();
        assert_eq!(pent, 12, "pentagon count != 12 for F={f}");
        let hex_ok = g.n1.iter().filter(|n| n.len() == 6).count();
        assert_eq!(pent + hex_ok, g.cells, "degree counts mismatch for F={f}");

        // Neighbor lists sorted (monotonic non-decreasing)
        for (i, n) in g.n1.iter().enumerate() {
            let mut prev = 0u32;
            for (k, &v) in n.iter().enumerate() {
                if k > 0 {
                    assert!(prev <= v, "n1 not sorted at cell {i}");
                }
                prev = v;
            }
        }
        for (i, n) in g.n2.iter().enumerate() {
            let mut prev = 0u32;
            for (k, &v) in n.iter().enumerate() {
                if k > 0 {
                    assert!(prev <= v, "n2 not sorted at cell {i}");
                }
                prev = v;
            }
        }
    }
}

#[test]
#[ignore]
fn area_partition_f64_diag() {
    let f = 64u32;
    let g = Grid::new(f);
    let sum_area: f64 = g.area.iter().map(|a| *a as f64).sum();
    let sphere = 4.0 * std::f64::consts::PI;
    let total_rel = rel_err(sum_area, sphere);
    println!("F=64 total area rel err: {:.3e}", total_rel);
    assert!(total_rel < 1e-6);
}

#[test]
fn cache_round_trip_minimal() -> Result<(), GridError> {
    let g = Grid::new(0);
    let ts = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_millis();
    let path = env::temp_dir().join(format!("grid_{}_{}.cache", std::process::id(), ts));
    g.save_cache(&path)?;
    let g2 = Grid::load_cache(&path)?;
    assert_eq!(g.cells, g2.cells);
    assert_eq!(g.pos_xyz, g2.pos_xyz);
    assert_eq!(g.latlon, g2.latlon);
    assert_eq!(g.area, g2.area);
    assert_eq!(g.n1, g2.n1);
    assert_eq!(g.n2, g2.n2);
    let _ = fs::remove_file(&path);
    Ok(())
}
