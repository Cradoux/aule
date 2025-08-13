//! Plate seeds, Euler poles, and per-cell velocities (T-030).

use crate::grid::Grid;
use rand::{RngCore, SeedableRng};

const RADIUS_M: f64 = 6_371_000.0;

/// Plate model results computed from a deterministic seed.
pub struct Plates {
    /// Plate seed unit vectors (f64)
    pub seeds: Vec<[f64; 3]>,
    /// Plate id per cell (u16)
    pub plate_id: Vec<u16>,
    /// Euler pole axis per plate (unit xyz)
    pub pole_axis: Vec<[f32; 3]>,
    /// Angular speed (rad/yr) per plate
    pub omega_rad_yr: Vec<f32>,
    /// Per-cell velocity components (east, north) in m/yr
    pub vel_en: Vec<[f32; 2]>,
}

impl Plates {
    /// Build plates and velocities for `num_plates` with a deterministic `seed`.
    pub fn new(grid: &Grid, num_plates: u32, seed: u64) -> Self {
        let seeds = farthest_point_seeds(grid, num_plates, seed ^ 0x0070_6c61_7465);
        let plate_id = assign_voronoi(grid, &seeds);
        let (pole_axis, omega_rad_yr) = euler_poles(seeds.len(), seed ^ 0xFACE_CAFE);
        let vel_en = velocities(grid, &plate_id, &pole_axis, &omega_rad_yr);
        Self { seeds, plate_id, pole_axis, omega_rad_yr, vel_en }
    }
}

/// Compute per-cell 3D surface velocities (m/yr) using current `plate_id` and plate kinematics.
pub fn velocity_field_m_per_yr(grid: &Grid, plates: &Plates, plate_id: &[u16]) -> Vec<[f32; 3]> {
    let mut vel = vec![[0.0f32, 0.0f32, 0.0f32]; grid.cells];
    for i in 0..grid.cells {
        let r = [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
        let pid = plate_id[i] as usize;
        let a = [
            plates.pole_axis[pid][0] as f64,
            plates.pole_axis[pid][1] as f64,
            plates.pole_axis[pid][2] as f64,
        ];
        let omega = plates.omega_rad_yr[pid] as f64;
        // linear velocity v = (omega * axis) × (R * r)
        let w = [a[0] * omega, a[1] * omega, a[2] * omega];
        let rp = [r[0] * RADIUS_M, r[1] * RADIUS_M, r[2] * RADIUS_M];
        let v = cross(w, rp);
        vel[i] = [v[0] as f32, v[1] as f32, v[2] as f32];
    }
    vel
}

/// Compute per-cell local east/north velocity components (m/yr) for the given `plate_id`.
pub fn velocity_en_m_per_yr(grid: &Grid, plates: &Plates, plate_id: &[u16]) -> Vec<[f32; 2]> {
    let mut vel_en = vec![[0.0f32, 0.0f32]; grid.cells];
    for i in 0..grid.cells {
        let p = [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
        let pid = plate_id[i] as usize;
        let a = [
            plates.pole_axis[pid][0] as f64,
            plates.pole_axis[pid][1] as f64,
            plates.pole_axis[pid][2] as f64,
        ];
        let omega = plates.omega_rad_yr[pid] as f64;
        let w = [a[0] * omega, a[1] * omega, a[2] * omega];
        let rp = [p[0] * RADIUS_M, p[1] * RADIUS_M, p[2] * RADIUS_M];
        let v = cross(w, rp);
        let east = normalize(cross([0.0, 0.0, 1.0], p));
        let north = normalize(cross(p, east));
        vel_en[i] = [dot(v, east) as f32, dot(v, north) as f32];
    }
    vel_en
}

/// Rotate a unit vector `r` about `axis` (unit) by angle `theta` (radians) using Rodrigues formula.
pub fn rotate_point(axis: [f32; 3], theta: f32, r: [f32; 3]) -> [f32; 3] {
    let k = [axis[0] as f64, axis[1] as f64, axis[2] as f64];
    let rr = [r[0] as f64, r[1] as f64, r[2] as f64];
    let ct = (theta as f64).cos();
    let st = (theta as f64).sin();
    let kxr = cross(k, rr);
    let kdotr = dot(k, rr);
    let term1 = [rr[0] * ct, rr[1] * ct, rr[2] * ct];
    let term2 = [kxr[0] * st, kxr[1] * st, kxr[2] * st];
    let term3 = [k[0] * kdotr * (1.0 - ct), k[1] * kdotr * (1.0 - ct), k[2] * kdotr * (1.0 - ct)];
    let out = [
        term1[0] + term2[0] + term3[0],
        term1[1] + term2[1] + term3[1],
        term1[2] + term2[2] + term3[2],
    ];
    // Normalize to unit length to mitigate drift
    let n = (out[0] * out[0] + out[1] * out[1] + out[2] * out[2]).sqrt();
    if n > 0.0 {
        [(out[0] / n) as f32, (out[1] / n) as f32, (out[2] / n) as f32]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn farthest_point_seeds(grid: &Grid, k: u32, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let first = (rng.next_u32() as usize) % grid.cells;

    let pos: Vec<[f64; 3]> =
        grid.pos_xyz.iter().map(|p| [p[0] as f64, p[1] as f64, p[2] as f64]).collect();

    let mut seeds_axes: Vec<[f64; 3]> = vec![pos[first]];
    let mut best_cos: Vec<f64> = vec![-1.0; grid.cells];
    for i in 0..grid.cells {
        best_cos[i] = dot(pos[i], pos[first]);
    }

    while (seeds_axes.len() as u32) < k {
        let mut min_cos = 1.0;
        let mut min_idx = 0usize;
        for (i, &c) in best_cos.iter().enumerate().take(grid.cells) {
            if c < min_cos {
                min_cos = c;
                min_idx = i;
            }
        }
        seeds_axes.push(pos[min_idx]);
        for i in 0..grid.cells {
            let c = dot(pos[i], pos[min_idx]);
            if c > best_cos[i] {
                best_cos[i] = c;
            }
        }
    }
    seeds_axes
}

fn assign_voronoi(grid: &Grid, seeds_axes: &[[f64; 3]]) -> Vec<u16> {
    const EPS_DOT: f64 = 1e-12;
    let pos: Vec<[f64; 3]> =
        grid.pos_xyz.iter().map(|p| [p[0] as f64, p[1] as f64, p[2] as f64]).collect();
    let mut plate_id = vec![0u16; grid.cells];
    for i in 0..grid.cells {
        let r = pos[i];
        let mut best_id: u16 = 0;
        let mut best_dot: f64 = -1.0;
        for (k, s) in seeds_axes.iter().enumerate() {
            let d = dot(r, *s).clamp(-1.0, 1.0);
            if d > best_dot + EPS_DOT || ((d - best_dot).abs() <= EPS_DOT && (k as u16) < best_id) {
                best_dot = d;
                best_id = k as u16;
            }
        }
        plate_id[i] = best_id;
    }
    plate_id
}

fn euler_poles(num_plates: usize, seed: u64) -> (Vec<[f32; 3]>, Vec<f32>) {
    // Deterministic poles: axis from seed index hash; omega from rng
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut axes: Vec<[f32; 3]> = Vec::with_capacity(num_plates);
    let mut omegas: Vec<f32> = Vec::with_capacity(num_plates);
    for i in 0..num_plates {
        // simple hash → axis on sphere
        let h = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let u = ((h ^ (h >> 33)) as f64 / u64::MAX as f64).clamp(0.0, 1.0);
        let v = ((h.rotate_left(13)) as f64 / u64::MAX as f64).clamp(0.0, 1.0);
        let theta = (2.0 * std::f64::consts::PI) * u;
        let z = 2.0 * v - 1.0;
        let r = (1.0 - z * z).sqrt();
        let axis = [r * theta.cos(), r * theta.sin(), z];
        axes.push([axis[0] as f32, axis[1] as f32, axis[2] as f32]);
        // omega rad/yr from rng in a small realistic range (e.g., up to ~1e-7 rad/yr)
        let omega = (rng.next_u32() as f64 / u32::MAX as f64) * 1.0e-7;
        omegas.push(omega as f32);
    }
    (axes, omegas)
}

fn velocities(
    grid: &Grid,
    plate_id: &[u16],
    pole_axis: &[[f32; 3]],
    omega_rad_yr: &[f32],
) -> Vec<[f32; 2]> {
    let mut vel = vec![[0.0f32, 0.0f32]; grid.cells];
    for i in 0..grid.cells {
        let p = [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
        let pid = plate_id[i] as usize;
        let a = [pole_axis[pid][0] as f64, pole_axis[pid][1] as f64, pole_axis[pid][2] as f64];
        let omega = omega_rad_yr[pid] as f64;
        // Angular velocity vector
        let w = [a[0] * omega, a[1] * omega, a[2] * omega];
        // Linear velocity v = w × (R * p)
        let rp = [p[0] * RADIUS_M, p[1] * RADIUS_M, p[2] * RADIUS_M];
        let v = cross(w, rp);
        // Project to local east/north basis at p
        let east = normalize(cross([0.0, 0.0, 1.0], p));
        let north = normalize(cross(p, east));
        let ve = dot(v, east) as f32;
        let vn = dot(v, north) as f32;
        vel[i] = [ve, vn];
    }
    vel
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n == 0.0 {
        return [0.0, 0.0, 0.0];
    }
    [v[0] / n, v[1] / n, v[2] / n]
}
