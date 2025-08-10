//! Plate seeds, Euler poles, and per-cell velocities (T-030).

use crate::grid::Grid;
use rand::{RngCore, SeedableRng};

const RADIUS_M: f64 = 6_371_000.0;

/// Plate model results.
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
    /// Build plates and velocities for N plates with a deterministic seed.
    pub fn new(grid: &Grid, num_plates: u32, seed: u64) -> Self {
        let seeds = farthest_point_seeds(grid, num_plates, seed ^ 0x0070_6c61_7465);
        let plate_id = assign_voronoi(grid, &seeds);
        let (pole_axis, omega_rad_yr) = euler_poles(seeds.len(), seed ^ 0xFACE_CAFE);
        let vel_en = velocities(grid, &plate_id, &pole_axis, &omega_rad_yr);
        Self { seeds, plate_id, pole_axis, omega_rad_yr, vel_en }
    }
}

fn farthest_point_seeds(grid: &Grid, k: u32, seed: u64) -> Vec<[f64; 3]> {
    // Deterministic RNG selection of the first seed
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let first = (rng.next_u32() as usize) % grid.cells;
    // unit vectors list of seeds
    let mut seeds_axes: Vec<[f64; 3]> = vec![[
        grid.pos_xyz[first][0] as f64,
        grid.pos_xyz[first][1] as f64,
        grid.pos_xyz[first][2] as f64,
    ]];
    // Precompute positions in f64
    let mut pos: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
    for p in &grid.pos_xyz {
        pos.push([p[0] as f64, p[1] as f64, p[2] as f64]);
    }
    let mut best_cos: Vec<f64> = vec![-1.0; grid.cells];
    // Initialize best_cos with first seed
    for i in 0..grid.cells {
        best_cos[i] = dot(pos[i], pos[first]);
    }
    while (seeds_axes.len() as u32) < k {
        // pick the cell with the smallest cosine (largest angle) to its nearest seed
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
