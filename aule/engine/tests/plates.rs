use engine::{grid::Grid, plates::Plates};

#[test]
fn determinism_and_voronoi_check() {
    let g = Grid::new(16);
    let p1 = Plates::new(&g, 8, 42);
    let p2 = Plates::new(&g, 8, 42);
    assert_eq!(p1.plate_id, p2.plate_id);
    assert_eq!(p1.pole_axis, p2.pole_axis);
    assert_eq!(p1.omega_rad_yr, p2.omega_rad_yr);
    assert_eq!(p1.vel_en, p2.vel_en);

    // Voronoi: each cell's assigned seed is indeed nearest by angle among seeds
    // Reconstruct seeds from plate ids by picking first index per plate
    let mut rep: Vec<usize> = vec![usize::MAX; 8];
    for (i, &pid) in p1.plate_id.iter().enumerate() {
        let k = pid as usize;
        if rep[k] == usize::MAX { rep[k] = i; }
    }
    for (i, &pid) in p1.plate_id.iter().enumerate() {
        let pos = [g.pos_xyz[i][0] as f64, g.pos_xyz[i][1] as f64, g.pos_xyz[i][2] as f64];
        let mut best = -1.0;
        let mut best_k = 0usize;
        for (k, &r) in rep.iter().enumerate() {
            let rpos = [g.pos_xyz[r][0] as f64, g.pos_xyz[r][1] as f64, g.pos_xyz[r][2] as f64];
            let c = pos[0]*rpos[0] + pos[1]*rpos[1] + pos[2]*rpos[2];
            if c > best { best = c; best_k = k; }
        }
        assert_eq!(best_k as u32, pid);
    }
}


