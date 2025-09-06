use engine::{boundaries::Boundaries, geo, grid::Grid};

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
    // edge_kin populated and aligned
    assert_eq!(b1.edge_kin.len(), b1.edges.len());

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

#[test]
fn edge_kin_consistency_spotcheck() {
    let f = 16u32;
    let g = Grid::new(f);
    let plates = engine::plates::Plates::new(&g, 8, 7);
    let b = Boundaries::classify(&g, &plates.plate_id, &plates.vel_en, 0.005);
    assert_eq!(b.edge_kin.len(), b.edges.len());
    for i in 0..b.edges.len() {
        let (u, v, class) = b.edges[i];
        let ek = &b.edge_kin[i];
        assert_eq!(ek.u, u);
        assert_eq!(ek.v, v);
        assert_eq!(ek.class as u8, class);
        // recompute dv·n_hat and dv·t_hat quickly
        let ru = [
            g.pos_xyz[u as usize][0] as f64,
            g.pos_xyz[u as usize][1] as f64,
            g.pos_xyz[u as usize][2] as f64,
        ];
        let rv = [
            g.pos_xyz[v as usize][0] as f64,
            g.pos_xyz[v as usize][1] as f64,
            g.pos_xyz[v as usize][2] as f64,
        ];
        let rm = geo::normalize([ru[0] + rv[0], ru[1] + rv[1], ru[2] + rv[2]]);
        let ggc = geo::normalize(geo::cross(ru, rv));
        let t_hat = geo::normalize(geo::cross(ggc, rm));
        let mut n_hat = geo::normalize(geo::cross(rm, t_hat));
        if geo::dot(n_hat, [rv[0] - ru[0], rv[1] - ru[1], rv[2] - ru[2]]) < 0.0 {
            n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
        }
        let en_to_w = |idx: usize| {
            let r = [g.pos_xyz[idx][0] as f64, g.pos_xyz[idx][1] as f64, g.pos_xyz[idx][2] as f64];
            let (e, n) = geo::local_basis(r);
            [
                (plates.vel_en[idx][0] as f64) * e[0] + (plates.vel_en[idx][1] as f64) * n[0],
                (plates.vel_en[idx][0] as f64) * e[1] + (plates.vel_en[idx][1] as f64) * n[1],
                (plates.vel_en[idx][0] as f64) * e[2] + (plates.vel_en[idx][1] as f64) * n[2],
            ]
        };
        let vu = en_to_w(u as usize);
        let vv = en_to_w(v as usize);
        let dv = [vu[0] - vv[0], vu[1] - vv[1], vu[2] - vv[2]];
        let n = geo::dot(dv, n_hat) as f32;
        let t = geo::dot(dv, t_hat) as f32;
        assert!((n - ek.n_m_per_yr).abs() < 1e-9);
        assert!((t - ek.t_m_per_yr).abs() < 1e-9);
    }
}
