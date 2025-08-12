use engine::{boundaries::Boundaries, grid::Grid, plates::Plates, transforms};

fn mean(xs: &[f32]) -> f32 {
    if xs.is_empty() {
        0.0
    } else {
        xs.iter().sum::<f32>() / xs.len() as f32
    }
}

#[test]
fn determinism_sign_and_idempotence() {
    let f: u32 = 16;
    let g = Grid::new(f);
    let plates = Plates::new(&g, 2, 2025);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);

    let mut depth0 = vec![3000.0_f32; g.cells];
    let p = transforms::TransformParams {
        tau_open_m_per_yr: 0.005,
        min_tangential_m_per_yr: 0.005,
        basin_half_width_km: 25.0,
        ridge_like_uplift_m: -200.0,
        basin_deepen_m: 400.0,
    };
    let (m0, s0) =
        transforms::apply_transforms(&g, &b, &plates.plate_id, &plates.vel_en, &mut depth0, p);

    // Determinism/idempotence
    let mut depth1 = vec![3000.0_f32; g.cells];
    let (m1, s1) =
        transforms::apply_transforms(&g, &b, &plates.plate_id, &plates.vel_en, &mut depth1, p);
    assert_eq!(s0.pull_apart_cells, s1.pull_apart_cells);
    assert_eq!(s0.restraining_cells, s1.restraining_cells);
    assert_eq!(m0.pull_apart, m1.pull_apart);
    assert_eq!(m0.restraining, m1.restraining);

    let mut depth2 = depth1.clone();
    let _ = transforms::apply_transforms(&g, &b, &plates.plate_id, &plates.vel_en, &mut depth2, p);
    assert_eq!(depth1, depth2);

    // Sign checks if non-empty
    if s0.pull_apart_cells > 0 {
        let mut deltas = Vec::new();
        for (i, &m) in m0.pull_apart.iter().enumerate() {
            if m {
                deltas.push(depth0[i] - 3000.0);
            }
        }
        assert!(mean(&deltas) > 0.0);
    }
    if s0.restraining_cells > 0 {
        let mut deltas = Vec::new();
        for (i, &m) in m0.restraining.iter().enumerate() {
            if m {
                deltas.push(depth0[i] - 3000.0);
            }
        }
        assert!(mean(&deltas) < 0.0);
    }
}
