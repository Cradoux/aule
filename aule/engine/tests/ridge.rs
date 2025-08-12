use engine::{boundaries::Boundaries, grid::Grid, plates::Plates, ridge};

#[test]
fn births_and_fringe_correctness_and_idempotence() {
    let f: u32 = 16;
    let g = Grid::new(f);
    // Deterministic plates: choose 2 plates for a simple ridge field
    let plates = Plates::new(&g, 2, 42);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);

    let mut ages = vec![10.0f32; g.cells];
    let p = ridge::RidgeParams { fringe_age_myr: 0.2 };
    let stats1 = ridge::apply_ridge(&g, &b, &mut ages, p);

    // All ridge cells must be age 0
    let mut ridge_count = 0usize;
    for &(u, v, class) in &b.edges {
        if class == 1 {
            if ages[u as usize] == 0.0 {
                ridge_count += 1;
            }
            if ages[v as usize] == 0.0 {
                ridge_count += 1;
            }
        }
    }
    // Each divergent edge contributes up to two ridge cells; counts may overlap, so >= births
    assert!(ridge_count as u32 >= stats1.births);

    // Fringe neighbors must be <= fringe cap only for neighbors of ridge cells
    let mut neighbor_mask = vec![false; g.cells];
    let mut is_ridge_cell = vec![false; g.cells];
    for &(u, v, class) in &b.edges {
        if class == 1 {
            is_ridge_cell[u as usize] = true;
            is_ridge_cell[v as usize] = true;
        }
    }
    for (u, &is_r) in is_ridge_cell.iter().enumerate() {
        if is_r {
            for &n in &g.n1[u] {
                neighbor_mask[n as usize] = true;
            }
        }
    }
    for i in 0..g.cells {
        if neighbor_mask[i] && ages[i] != 0.0 {
            assert!(ages[i] <= p.fringe_age_myr + 1e-6);
        }
    }

    // Idempotence: applying again should not change counts further
    let stats2 = ridge::apply_ridge(&g, &b, &mut ages, p);
    assert_eq!(stats2.births, 0);
    assert_eq!(stats2.fringe, 0);

    // Non-negativity and array length unchanged
    assert_eq!(ages.len(), g.cells);
    assert!(ages.iter().all(|&a| a >= 0.0));
}

#[test]
fn no_fringe_when_disabled() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 2, 7);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);

    let mut ages = vec![5.0f32; g.cells];
    let p = ridge::RidgeParams { fringe_age_myr: 0.0 };
    let stats = ridge::apply_ridge(&g, &b, &mut ages, p);
    assert!(stats.fringe == 0);
    // Births still occur
    assert!(stats.births > 0);
}
