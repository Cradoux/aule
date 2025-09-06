//! Plate seeds, Euler poles, and per-cell velocities (T-030).

use crate::grid::Grid;
use rand::{RngCore, SeedableRng};

const RADIUS_M: f64 = 6_371_000.0;

/// Sentinel used to denote an invalid/missing plate id.
pub const INVALID_PLATE_ID: u16 = u16::MAX;

/// Categorical plate type used to gate boundary physics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlateKind {
    /// Plate whose interior area has a majority of high-`C` cells (continent-dominated)
    Continental,
    /// Plate whose interior area is mostly low-`C` cells (ocean-dominated)
    Oceanic,
}

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
    /// Per-plate kind derived from initial continental mask
    pub kind: Vec<PlateKind>,
}

impl Plates {
    /// Build plates and velocities for `num_plates` with a deterministic `seed`.
    pub fn new(grid: &Grid, num_plates: u32, seed: u64) -> Self {
        let seeds = farthest_point_seeds(grid, num_plates, seed ^ 0x0070_6c61_7465);
        let plate_id = assign_voronoi(grid, &seeds);
        let (pole_axis, mut omega_rad_yr) = euler_poles(seeds.len(), seed ^ 0xFACE_CAFE);
        // Clamp angular speeds to a reasonable max linear speed (e.g., 100 mm/yr)
        let vmax_m_per_yr: f64 = 100.0e-3; // 100 mm/yr in m/yr
        let omega_max = (vmax_m_per_yr / RADIUS_M) as f32; // rad/yr
        for w in &mut omega_rad_yr {
            let s = (*w).abs().min(omega_max);
            *w = (*w).signum() * s;
        }
        let vel_en = velocities(grid, &plate_id, &pole_axis, &omega_rad_yr);
        Self { seeds, plate_id, pole_axis, omega_rad_yr, vel_en, kind: Vec::new() }
    }

    /// Advance rigid plate mosaics by rotating Voronoi seeds around each plate's Euler pole.
    /// Seeds are rotated by theta = omega * dt_years (radians). After rotation, rebuild
    /// `plate_id` and `vel_en` to keep the mosaic and velocities consistent.
    pub fn advance_rigid(&mut self, grid: &Grid, dt_years: f64) {
        let nplates = self.pole_axis.len().min(self.omega_rad_yr.len());
        for p in 0..nplates {
            // Convert axis to f64
            let axis = [
                self.pole_axis[p][0] as f64,
                self.pole_axis[p][1] as f64,
                self.pole_axis[p][2] as f64,
            ];
            let theta = (self.omega_rad_yr[p] as f64) * dt_years;
            // Rotate each seed assigned to this plate: seeds are one per plate by construction
            if p < self.seeds.len() {
                self.seeds[p] = rotate_about_axis_f64(self.seeds[p], axis, theta);
            }
        }
        // Rebuild mosaic and velocities
        self.plate_id = assign_voronoi(grid, &self.seeds);
        self.vel_en = velocities(grid, &self.plate_id, &self.pole_axis, &self.omega_rad_yr);
    }

    /// Maximum absolute angular speed (rad/yr) among plates.
    pub fn max_abs_omega_rad_yr(&self) -> f64 {
        let mut m = 0.0_f64;
        for &w in &self.omega_rad_yr {
            let a = (w as f64).abs();
            if a > m {
                m = a;
            }
        }
        m
    }

    /// Spawn a new plate near `cell` by adding a seed/pole and reassigning a small local patch.
    ///
    /// Deterministic given (grid positions, current plate count, seed). Keeps existing plate indices stable.
    /// Returns the new plate id.
    pub fn spawn_plate_at(&mut self, grid: &Grid, cell: usize, seed: u64, ring_r: u32) -> u16 {
        let new_pid = self.seeds.len() as u16;
        // Seed axis from target cell position (normalized)
        let r = [
            grid.pos_xyz[cell][0] as f64,
            grid.pos_xyz[cell][1] as f64,
            grid.pos_xyz[cell][2] as f64,
        ];
        let r_n = normalize(r);
        self.seeds.push(r_n);
        // Simple deterministic pole axis and omega derived from seed and index
        let h = seed ^ (new_pid as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let axis = hash_to_unit(h);
        self.pole_axis.push([axis[0] as f32, axis[1] as f32, axis[2] as f32]);
        // Small angular speed in a plausible range, deterministic from hash bits
        let omega = ((h.rotate_left(13) as f64) / (u64::MAX as f64)) * 1.0e-7;
        self.omega_rad_yr.push(omega as f32);
        // Reassign a local patch to the new plate id using BFS up to ring_r rings
        let mut visited: Vec<bool> = vec![false; grid.cells];
        let mut q: std::collections::VecDeque<(usize, u32)> = std::collections::VecDeque::new();
        q.push_back((cell.min(grid.cells - 1), 0));
        visited[cell.min(grid.cells - 1)] = true;
        while let Some((u, d)) = q.pop_front() {
            self.plate_id[u] = new_pid;
            if d >= ring_r {
                continue;
            }
            for &vn in &grid.n1[u] {
                let v = vn as usize;
                if !visited[v] {
                    visited[v] = true;
                    q.push_back((v, d + 1));
                }
            }
        }
        // Refresh velocities for fallback path
        self.vel_en = velocities(grid, &self.plate_id, &self.pole_axis, &self.omega_rad_yr);
        new_pid
    }

    /// Soft-retire a plate by reassigning its cells to `merge_into` and zeroing its omega.
    /// Indices remain stable; velocities are rebuilt. No-op if ids are equal or out of range.
    pub fn retire_plate_soft(&mut self, grid: &Grid, plate_id_to_retire: u16, merge_into: u16) {
        let nplates = self.pole_axis.len() as u16;
        if plate_id_to_retire >= nplates
            || merge_into >= nplates
            || plate_id_to_retire == merge_into
        {
            return;
        }
        for pid in &mut self.plate_id {
            if *pid == plate_id_to_retire {
                *pid = merge_into;
            }
        }
        // Zero rotation so retired plate contributes no motion even if referenced
        self.omega_rad_yr[plate_id_to_retire as usize] = 0.0;
        self.vel_en = velocities(grid, &self.plate_id, &self.pole_axis, &self.omega_rad_yr);
    }
}

/// Compute per-cell 3D surface velocities (m/yr) using current `plate_id` and plate kinematics.
pub fn velocity_field_m_per_yr(grid: &Grid, plates: &Plates, plate_id: &[u16]) -> Vec<[f32; 3]> {
    let n = grid.cells;
    // Precompute per-plate angular velocity vectors w = omega * axis (f64)
    let wvec: Vec<[f64; 3]> = plates
        .pole_axis
        .iter()
        .zip(plates.omega_rad_yr.iter())
        .map(|(a, &w)| {
            let ww = w as f64;
            [a[0] as f64 * ww, a[1] as f64 * ww, a[2] as f64 * ww]
        })
        .collect();
    let threads = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).clamp(1, 16);
    let chunk = ((n + threads - 1) / threads).max(1);
    let mut out = vec![[0.0f32; 3]; n];
    std::thread::scope(|scope| {
        for t in 0..threads {
            let start = t * chunk;
            if start >= n {
                break;
            }
            let end = (start + chunk).min(n);
            let pos = grid.pos_xyz[start..end].to_vec();
            let pid_slice = plate_id[start..end].to_vec();
            let wvec = wvec.clone();
            let handle = scope.spawn(move || {
                let mut loc = vec![[0.0f32; 3]; pos.len()];
                for k in 0..pos.len() {
                    let r = [pos[k][0] as f64, pos[k][1] as f64, pos[k][2] as f64];
                    let pid = pid_slice[k] as usize;
                    let w = wvec[pid];
                    let rp = [r[0] * RADIUS_M, r[1] * RADIUS_M, r[2] * RADIUS_M];
                    let v = cross(w, rp);
                    loc[k] = [v[0] as f32, v[1] as f32, v[2] as f32];
                }
                (start, loc)
            });
            let (s0, loc) = handle.join().unwrap();
            out[s0..s0 + loc.len()].copy_from_slice(&loc);
        }
    });
    out
}

/// Compute per-cell local east/north velocity components (m/yr) for the given `plate_id`.
pub fn velocity_en_m_per_yr(grid: &Grid, plates: &Plates, plate_id: &[u16]) -> Vec<[f32; 2]> {
    let n = grid.cells;
    // Precompute per-plate angular velocity vectors w = omega * axis (f64)
    let wvec: Vec<[f64; 3]> = plates
        .pole_axis
        .iter()
        .zip(plates.omega_rad_yr.iter())
        .map(|(a, &w)| {
            let ww = w as f64;
            [a[0] as f64 * ww, a[1] as f64 * ww, a[2] as f64 * ww]
        })
        .collect();
    let threads = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).clamp(1, 16);
    let chunk = ((n + threads - 1) / threads).max(1);
    let mut out = vec![[0.0f32; 2]; n];
    std::thread::scope(|scope| {
        for t in 0..threads {
            let start = t * chunk;
            if start >= n {
                break;
            }
            let end = (start + chunk).min(n);
            let pos = grid.pos_xyz[start..end].to_vec();
            let east = grid.east_hat[start..end].to_vec();
            let north = grid.north_hat[start..end].to_vec();
            let pid_slice = plate_id[start..end].to_vec();
            let wvec = wvec.clone();
            let handle = scope.spawn(move || {
                let mut loc = vec![[0.0f32; 2]; pos.len()];
                for k in 0..pos.len() {
                    let p = [pos[k][0] as f64, pos[k][1] as f64, pos[k][2] as f64];
                    let pid = pid_slice[k] as usize;
                    let w = wvec[pid];
                    let rp = [p[0] * RADIUS_M, p[1] * RADIUS_M, p[2] * RADIUS_M];
                    let v = cross(w, rp);
                    let eh = [east[k][0] as f64, east[k][1] as f64, east[k][2] as f64];
                    let nh = [north[k][0] as f64, north[k][1] as f64, north[k][2] as f64];
                    loc[k] = [dot(v, eh) as f32, dot(v, nh) as f32];
                }
                (start, loc)
            });
            let (s0, loc) = handle.join().unwrap();
            out[s0..s0 + loc.len()].copy_from_slice(&loc);
        }
    });
    out
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

/// Heal invalid plate ids in-place by assigning each invalid cell to the nearest valid plate.
///
/// Strategy: multi-source BFS starting from all valid cells; propagate their plate ids into
/// neighbouring invalid regions using the 1-ring adjacency until all cells are labeled.
pub fn heal_plate_ids(grid: &Grid, plate_id: &mut [u16]) {
    let n = grid.cells.min(plate_id.len());
    if n == 0 {
        return;
    }
    // Collect valid seeds
    let mut dist: Vec<i32> = vec![-1; n];
    let mut q: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    for i in 0..n {
        if plate_id[i] != INVALID_PLATE_ID {
            dist[i] = 0;
            q.push_back(i);
        }
    }
    if q.is_empty() {
        // No valid seeds: assign a default plate id 0 to all cells
        for pid in plate_id.iter_mut().take(n) {
            *pid = 0;
        }
        return;
    }
    while let Some(u) = q.pop_front() {
        let pid_u = plate_id[u];
        for &vn in &grid.n1[u] {
            let v = vn as usize;
            if v >= n {
                continue;
            }
            if plate_id[v] == INVALID_PLATE_ID && dist[v] < 0 {
                plate_id[v] = pid_u;
                dist[v] = dist[u] + 1;
                q.push_back(v);
            }
        }
    }
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

#[inline]
fn hash_to_unit(h: u64) -> [f64; 3] {
    let u = ((h ^ (h >> 33)) as f64 / u64::MAX as f64).clamp(0.0, 1.0);
    let v = ((h.rotate_left(13)) as f64 / u64::MAX as f64).clamp(0.0, 1.0);
    let theta = (2.0 * std::f64::consts::PI) * u;
    let z = 2.0 * v - 1.0;
    let r = (1.0 - z * z).sqrt();
    [r * theta.cos(), r * theta.sin(), z]
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
        // Project to local east/north basis at p via geo::local_basis
        let (east, north) = crate::geo::local_basis(p);
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

/// Derive per-plate kinds from area-weighted continental fraction.
pub fn derive_kinds(
    grid: &Grid,
    plate_id: &[u16],
    c: &[f32],
    frac_thresh: f32,
    c_thresh: f32,
) -> Vec<PlateKind> {
    let nplates = plate_id.iter().copied().map(|p| p as usize).max().unwrap_or(0) + 1;
    let mut a_tot = vec![0.0f64; nplates];
    let mut a_con = vec![0.0f64; nplates];
    for i in 0..grid.cells.min(plate_id.len()).min(c.len()) {
        let pid = plate_id[i] as usize;
        let a = grid.area[i] as f64;
        a_tot[pid] += a;
        if c[i] > c_thresh {
            a_con[pid] += a;
        }
    }
    let mut kinds = vec![PlateKind::Oceanic; nplates];
    for pid in 0..nplates {
        let frac = if a_tot[pid] > 0.0 { a_con[pid] / a_tot[pid] } else { 0.0 } as f32;
        kinds[pid] = if frac > frac_thresh { PlateKind::Continental } else { PlateKind::Oceanic };
    }
    kinds
}

/// Rodrigues rotation of a unit vector `u` about unit `axis` by `theta` radians (f64).
pub(crate) fn rotate_about_axis_f64(u: [f64; 3], axis: [f64; 3], theta: f64) -> [f64; 3] {
    // Normalize axis defensively
    let a = normalize(axis);
    let (s, c) = theta.sin_cos();
    // a*(a·u)*(1-c) + u*c + (a×u)*s
    let adotu = a[0] * u[0] + a[1] * u[1] + a[2] * u[2];
    let axu = cross(a, u);
    let term1 = [a[0] * adotu * (1.0 - c), a[1] * adotu * (1.0 - c), a[2] * adotu * (1.0 - c)];
    let term2 = [u[0] * c, u[1] * c, u[2] * c];
    let term3 = [axu[0] * s, axu[1] * s, axu[2] * s];
    let out = [
        term1[0] + term2[0] + term3[0],
        term1[1] + term2[1] + term3[1],
        term1[2] + term2[2] + term3[2],
    ];
    normalize(out)
}
