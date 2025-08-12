use engine::{boundaries::Boundaries, grid::Grid};

#[test]
fn determinism_symmetry_basic() {
    let f = 16u32;
    let g = Grid::new(f);
    // Build a deterministic plate/velocity field using Plates module
    let plates = engine::plates::Plates::new(&g, 8, 42);
    let b1 = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    let b2 = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    assert_eq!(b1.b, b2.b);
    assert_eq!(b1.edges, b2.edges);
    assert_eq!(b1.stats, b2.stats);

    // Symmetry: each edge sets bits on both endpoints
    for (u, v, class) in b1.edges {
        let mask = match class {
            1 => 1,
            2 => 2,
            3 => 4,
            _ => 0,
        };
        assert_ne!(mask, 0);
        assert!(b1.b[u as usize] & mask != 0);
        assert!(b1.b[v as usize] & mask != 0);
    }
}

#[test]
fn threshold_monotonicity() {
    let f = 16u32;
    let g = Grid::new(f);
    let plates = engine::plates::Plates::new(&g, 8, 123);
    let hi = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.02);
    let lo = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.002);
    // Lower tau should not flip div<->conv for the same edge; it can only convert some transforms
    // into div or conv. Here we only assert counts are sensible and non-zero across classes.
    assert!(hi.stats.divergent + hi.stats.convergent + hi.stats.transform > 0);
    assert!(
        lo.stats.divergent + lo.stats.convergent + lo.stats.transform
            >= hi.stats.divergent + hi.stats.convergent
    );
}
