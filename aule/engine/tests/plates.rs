use engine::{grid::Grid, plates::{Plates, INVALID_PLATE_ID, heal_plate_ids}};

#[test]
fn determinism_and_voronoi_check() {
    let g = Grid::new(16);
    let p1 = Plates::new(&g, 8, 42);
    let p2 = Plates::new(&g, 8, 42);
    assert_eq!(p1.plate_id, p2.plate_id);
    assert_eq!(p1.pole_axis, p2.pole_axis);
    assert_eq!(p1.omega_rad_yr, p2.omega_rad_yr);
    assert_eq!(p1.vel_en, p2.vel_en);

    // Voronoi: each cell's assigned seed is nearest by max dot with tie-breaker
    let rep_axes: &[[f64; 3]] = &p1.seeds;
    const EPS_DOT: f64 = 1e-12;
    for (i, &pid) in p1.plate_id.iter().enumerate() {
        let pos = [g.pos_xyz[i][0] as f64, g.pos_xyz[i][1] as f64, g.pos_xyz[i][2] as f64];
        let mut best_dot = -1.0f64;
        let mut best_id: u16 = 0;
        for (k, s) in rep_axes.iter().enumerate() {
            let c = (pos[0] * s[0] + pos[1] * s[1] + pos[2] * s[2]).clamp(-1.0, 1.0);
            if c > best_dot + EPS_DOT || ((c - best_dot).abs() <= EPS_DOT && (k as u16) < best_id) {
                best_dot = c;
                best_id = k as u16;
            }
        }
        assert_eq!(best_id, pid);
    }
}

#[test]
fn heal_plate_ids_fills_invalid_regions() {
    let g = Grid::new(8);
    let mut p = Plates::new(&g, 4, 7);
    // Invalidate a contiguous patch (seed center and its 1-ring)
    let center = g.cells / 2;
    p.plate_id[center] = INVALID_PLATE_ID;
    for &n in &g.n1[center] {
        p.plate_id[n as usize] = INVALID_PLATE_ID;
    }
    // Ensure we have at least one valid elsewhere
    assert!(p.plate_id.iter().any(|&x| x != INVALID_PLATE_ID));
    // Heal
    heal_plate_ids(&g, &mut p.plate_id);
    // No invalids remain
    assert!(p.plate_id.iter().all(|&x| x != INVALID_PLATE_ID));
}
