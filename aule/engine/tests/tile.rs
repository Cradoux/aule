use engine::grid::tile::Tiling;
use engine::grid::Grid;

#[test]
fn tiler_coverage_and_halo_properties() {
    let g = Grid::new(16);
    let t = Tiling::new(&g, 2, 8192);

    // Coverage & exclusivity: every cell owned by exactly one interior
    assert_eq!(t.owner.len(), g.cells);
    let mut count_owned = 0usize;
    for (i, &own) in t.owner.iter().enumerate() {
        assert!(own != u32::MAX, "cell {i} has no owner");
        count_owned += 1;
    }
    assert_eq!(count_owned, g.cells);

    // Interiors disjoint
    let mut seen = vec![false; g.cells];
    for tile in &t.tiles {
        for &u in &tile.interior {
            let ui = u as usize;
            assert!(!seen[ui], "cell {ui} appears in two interiors");
            seen[ui] = true;
        }
    }

    // Halo correctness: for each tile and interior cell, 2-ring subset is in interior ∪ halo
    for tile in &t.tiles {
        let interior_set: std::collections::BTreeSet<u32> = tile.interior.iter().copied().collect();
        let halo_set: std::collections::BTreeSet<u32> = tile.halo.iter().copied().collect();
        for &u in &tile.interior {
            // 1-ring
            for &v in &g.n1[u as usize] {
                assert!(interior_set.contains(&v) || halo_set.contains(&v));
            }
            // 2-ring via neighbors of neighbors
            for &v in &g.n1[u as usize] {
                for &w in &g.n1[v as usize] {
                    if w != u && !g.n1[u as usize].contains(&w) {
                        assert!(interior_set.contains(&w) || halo_set.contains(&w));
                    }
                }
            }
        }
    }

    // Adjacency symmetric and unique
    for (i, a) in t.tiles.iter().enumerate() {
        for &b in &a.neighbors {
            let btile = &t.tiles[b as usize];
            assert!(btile.neighbors.binary_search(&(i as u32)).is_ok());
        }
        // uniqueness is implied by sorted set build
    }
}

#[test]
fn tiler_quick_f32() {
    let g = Grid::new(32);
    let t = Tiling::new(&g, 2, 8192);
    // Sanity: with per-face blocking and target >> per-face cells, expect one block per face (20)
    let f = 32usize;
    let per_face = (f + 1) * (f + 2) / 2; // 561
    assert!(8192 >= per_face);
    assert_eq!(t.tiles.len(), 20);
}
