//! Semi-Lagrangian backtrace helpers for categorical and scalar fields (MVP).

use crate::grid::Grid;

/// Backtrace from unit vector `r` by `dt_myr` along velocity `v_m_per_yr` (3D) on the sphere.
/// Returns a unit vector approximating the source position.
pub fn backtrace_r(r: [f32; 3], v_m_per_yr: [f32; 3], dt_myr: f64) -> [f32; 3] {
    // Convert displacement to an angular small rotation about axis perpendicular to r and v
    // For MVP, use simple Euler: r_src ≈ normalize(r - (dt * v / R))
    const R_EARTH_M: f64 = 6_371_000.0;
    let dt_yr = (dt_myr * 1.0e6) as f32;
    let dr = [
        -v_m_per_yr[0] * dt_yr / (R_EARTH_M as f32),
        -v_m_per_yr[1] * dt_yr / (R_EARTH_M as f32),
        -v_m_per_yr[2] * dt_yr / (R_EARTH_M as f32),
    ];
    let rr = [r[0] + dr[0], r[1] + dr[1], r[2] + dr[2]];
    let n = (rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]).sqrt();
    if n > 0.0 {
        [rr[0] / n, rr[1] / n, rr[2] / n]
    } else {
        r
    }
}

/// Find nearest neighbor cell to unit vector `r_src` using 1-ring + 2-ring of the provided `start_i`.
fn nearest_idx_local(grid: &Grid, start_i: usize, r_src: [f32; 3]) -> usize {
    let mut best_i = start_i;
    let mut best_d2 = f32::INFINITY;
    let check = |i: usize, best_i: &mut usize, best_d2: &mut f32| {
        let dx = grid.pos_xyz[i][0] - r_src[0];
        let dy = grid.pos_xyz[i][1] - r_src[1];
        let dz = grid.pos_xyz[i][2] - r_src[2];
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < *best_d2 || (d2 - *best_d2).abs() <= f32::EPSILON && i < *best_i {
            *best_d2 = d2;
            *best_i = i;
        }
    };
    check(start_i, &mut best_i, &mut best_d2);
    for &n in &grid.n1[start_i] {
        check(n as usize, &mut best_i, &mut best_d2);
    }
    for &n in &grid.n2[start_i] {
        check(n as usize, &mut best_i, &mut best_d2);
    }
    best_i
}

/// Advect categorical `plate_id` via nearest-neighbor semi-Lagrangian backtrace.
pub fn advect_plate_id(
    grid: &Grid,
    vel_m_per_yr: &[[f32; 3]],
    dt_myr: f64,
    plate_id_in: &[u16],
    plate_id_out: &mut [u16],
) {
    let n = grid.cells.min(vel_m_per_yr.len()).min(plate_id_in.len()).min(plate_id_out.len());
    for i in 0..n {
        let r = grid.pos_xyz[i];
        let v = vel_m_per_yr[i];
        let r_src = backtrace_r(r, v, dt_myr);
        // 1) Fast local guess using start cell's 1-2 ring
        let mut j = nearest_idx_local(grid, i, r_src);
        // 2) Robustify via bounded hill-climb over neighbour rings to reduce missed matches
        let mut best_d2 = {
            let dx = grid.pos_xyz[j][0] - r_src[0];
            let dy = grid.pos_xyz[j][1] - r_src[1];
            let dz = grid.pos_xyz[j][2] - r_src[2];
            dx * dx + dy * dy + dz * dz
        };
        // Up to 12 iterations of greedy improvement
        for _ in 0..24 {
            let mut improved = false;
            // Check 1-ring and 2-ring around current best
            let try_idx =
                |k: usize, j_cur: &mut usize, d2_best: &mut f32, improved_flag: &mut bool| {
                    let p = grid.pos_xyz[k];
                    let dx = p[0] - r_src[0];
                    let dy = p[1] - r_src[1];
                    let dz = p[2] - r_src[2];
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 < *d2_best || ((d2 - *d2_best).abs() <= f32::EPSILON && k < *j_cur) {
                        *d2_best = d2;
                        *j_cur = k;
                        *improved_flag = true;
                    }
                };
            for &k in &grid.n1[j] {
                try_idx(k as usize, &mut j, &mut best_d2, &mut improved);
            }
            for &k in &grid.n2[j] {
                try_idx(k as usize, &mut j, &mut best_d2, &mut improved);
            }
            if !improved {
                break;
            }
        }
        plate_id_out[i] = plate_id_in[j];
    }
}
