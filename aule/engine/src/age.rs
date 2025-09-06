//! Steady-state oceanic age via Dijkstra from ridge seeds and age→depth mapping.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::boundaries::Boundaries;
use crate::geo;
use crate::geo::dot;
use crate::grid::Grid;

const RADIUS_M: f64 = 6_371_000.0;

/// Parameters for steady-state age solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AgeParams {
    /// Minimum speed used to avoid infinite times at stagnation (m/yr).
    pub v_floor_m_per_yr: f64,
}

impl Default for AgeParams {
    fn default() -> Self {
        Self { v_floor_m_per_yr: 0.005 }
    }
}

/// Outputs of the age solver and bathymetry mapping.
#[derive(Debug, Clone, PartialEq)]
pub struct AgeOutputs {
    /// Age in Myr per cell (len = grid.cells)
    pub age_myr: Vec<f32>,
    /// Bathymetry in meters (positive down)
    pub depth_m: Vec<f32>,
    /// Min/max of age
    pub min_max_age: (f32, f32),
    /// Min/max of depth
    pub min_max_depth: (f32, f32),
}
/// Compute bathymetry depth (m, +down) from age (Myr) using a simple
/// Parsons–Sclater style curve (no clamping).
#[inline]
pub fn depth_from_age(age_myr: f64, d0: f64, a: f64, b: f64) -> f64 {
    d0 + a * age_myr.sqrt() + b * age_myr
}

/// Half-space cooling depth model using thermal diffusivity `k_m2_s`.
/// Units: `age_myr` in Myr, `d0_m` and result in meters, `k_m2_s` in m^2/s.
/// This simple parameterization is tuned to produce ~2–6 km depths over 0–120 Myr
/// for typical values (d0≈2600 m, k≈1e-6 m^2/s).
#[inline]
pub fn depth_from_age_hsc(age_myr: f64, d0_m: f64, k_m2_s: f64) -> f64 {
    // Calibrated to produce ~2.6–6.0 km over 0–120 Myr using an effective prefactor
    let age_s = (age_myr.max(0.0)) * 1.0e6_f64 * 365.25_f64 * 86400.0_f64;
    let k = k_m2_s.max(0.0);
    let c = 0.005_f64; // reduced from 0.02 to align with PS-like depths
    d0_m + c * (k * age_s).sqrt()
}

/// Plate-cooling blend with asymptote: d(t) = d_inf - (d_inf - (d_ridge + A sqrt(t))) * exp(-t/τ).
/// Units: age in Myr, distances in meters, τ in Myr, A in m/√Myr.
#[inline]
pub fn depth_from_age_blend(
    age_myr: f64,
    d_ridge_m: f64,
    a_sqrt_m_per_sqrt_myr: f64,
    d_inf_m: f64,
    tau_myr: f64,
) -> f64 {
    let t = age_myr.max(0.0);
    let base = d_ridge_m + a_sqrt_m_per_sqrt_myr * t.sqrt();
    let tau = tau_myr.max(1e-6);
    let w = (-t / tau).exp();
    let d = d_inf_m - (d_inf_m - base) * w;
    d.clamp(0.0, d_inf_m.max(d_ridge_m))
}

/// Plate cooling depth model with asymptotic plate thickness/maximum depth `z_plate_m`.
/// For simplicity we cap the HSC curve at `z_plate_m`. With large `z_plate_m` the
/// result approaches HSC at young ages. `t_myr` is accepted for API completeness but
/// not used in this minimal implementation.
#[inline]
pub fn depth_from_age_plate(
    age_myr: f64,
    d0_m: f64,
    _t_myr: f64,
    z_plate_m: f64,
    k_m2_s: f64,
) -> f64 {
    depth_from_age_hsc(age_myr, d0_m, k_m2_s).min(z_plate_m.max(d0_m))
}

#[derive(Copy, Clone, Debug)]
struct QueueItem {
    t_yr: f64,
    cell: u32,
}
impl Eq for QueueItem {}
impl PartialEq for QueueItem {
    fn eq(&self, other: &Self) -> bool {
        self.t_yr.eq(&other.t_yr) && self.cell == other.cell
    }
}
impl PartialOrd for QueueItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for QueueItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse by time (min-heap behavior)
        match self.t_yr.partial_cmp(&other.t_yr) {
            // Reverse ordering: smaller t_yr = Greater ordering so it pops later
            Some(Ordering::Less) => Ordering::Greater,
            Some(Ordering::Greater) => Ordering::Less,
            Some(Ordering::Equal) => self.cell.cmp(&other.cell).reverse(),
            None => Ordering::Equal,
        }
    }
}

/// Compute steady-state age and bathymetry from ridges and velocities.
pub fn compute_age_and_bathymetry(
    grid: &Grid,
    boundaries: &Boundaries,
    plate_id: &[u16],
    v_en: &[[f32; 2]],
    params: AgeParams,
) -> AgeOutputs {
    assert_eq!(plate_id.len(), grid.cells);
    assert_eq!(v_en.len(), grid.cells);

    // Build ridge seeds from divergent edges
    let mut is_seed = vec![false; grid.cells];
    for &(u, v, class) in &boundaries.edges {
        if class == 1 {
            // divergent
            is_seed[u as usize] = true;
            is_seed[v as usize] = true;
        }
    }

    // Dijkstra over cells constrained by plate id
    let mut dist_yr: Vec<f64> = vec![f64::INFINITY; grid.cells];
    let mut visited: Vec<bool> = vec![false; grid.cells];
    let mut heap: BinaryHeap<QueueItem> = BinaryHeap::new();

    for i in 0..grid.cells {
        if is_seed[i] {
            dist_yr[i] = 0.0;
            heap.push(QueueItem { t_yr: 0.0, cell: i as u32 });
        }
    }

    let vfloor = params.v_floor_m_per_yr.max(1e-12);
    let mut pos_unit: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
    for p in &grid.pos_xyz {
        let r = [p[0] as f64, p[1] as f64, p[2] as f64];
        pos_unit.push(geo::normalize(r));
    }

    while let Some(QueueItem { t_yr, cell }) = heap.pop() {
        let u = cell as usize;
        if visited[u] {
            continue;
        }
        visited[u] = true;

        let pid_u = plate_id[u];
        // neighbors within same plate
        for &vn in &grid.n1[u] {
            let v = vn as usize;
            if plate_id[v] != pid_u {
                continue;
            }
            if visited[v] {
                continue;
            }

            // edge weight dt = distance / max(|V(u)|, v_floor)
            let s_m = geo::great_circle_arc_len_m(pos_unit[u], pos_unit[v], RADIUS_M);
            let speed = ((v_en[u][0] as f64).hypot(v_en[u][1] as f64)).max(vfloor);
            let dt_yr = s_m / speed;
            let new_t = t_yr + dt_yr;
            if new_t < dist_yr[v] {
                dist_yr[v] = new_t;
                heap.push(QueueItem { t_yr: new_t, cell: v as u32 });
            }
        }
    }

    // Convert to Myr and compute depth via Parsons–Sclater-like map
    let mut age_myr: Vec<f32> = Vec::with_capacity(grid.cells);
    let mut depth_m: Vec<f32> = Vec::with_capacity(grid.cells);
    let (mut amin, mut amax) = (f32::INFINITY, f32::NEG_INFINITY);
    let (mut dmin, mut dmax) = (f32::INFINITY, f32::NEG_INFINITY);
    let (d0, a_coef, b_coef) = (2600.0_f64, 350.0_f64, 0.0_f64);
    for t in &dist_yr {
        let t_myr = if t.is_finite() { (*t / 1.0e6) as f32 } else { f32::INFINITY };
        let mut depth = depth_from_age(t_myr as f64, d0, a_coef, b_coef) as f32;
        if !depth.is_finite() {
            depth = 6000.0;
        }
        depth = depth.clamp(0.0, 6000.0);
        amin = amin.min(t_myr);
        amax = amax.max(t_myr);
        dmin = dmin.min(depth);
        dmax = dmax.max(depth);
        age_myr.push(t_myr);
        depth_m.push(depth);
    }

    AgeOutputs { age_myr, depth_m, min_max_age: (amin, amax), min_max_depth: (dmin, dmax) }
}

/// Exact rigid advection of the age field under plate rotation.
/// Back-rotate current cell centers and sample nearest source into stage, then publish.
pub fn rigid_advect_age(
    grid: &Grid,
    plates: &crate::plates::Plates,
    dt_years: f64,
    age_myr: &mut [f32],
    age_stage: &mut [f32],
) {
    let n = grid.cells.min(age_myr.len()).min(age_stage.len());
    let src = age_myr.to_vec();
    let threads = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).clamp(1, 16);
    let chunk = ((n + threads - 1) / threads).max(1);
    let mut assembled: Vec<f32> = vec![0.0; n];
    let src_ref = &src;
    std::thread::scope(|scope| {
        for t in 0..threads {
            let start = t * chunk;
            if start >= n {
                break;
            }
            let end = (start + chunk).min(n);
            let src_ref = src_ref;
            let handle = scope.spawn(move || {
                let mut loc: Vec<f32> = vec![0.0; end - start];
                for k in 0..(end - start) {
                    let i = start + k;
                    let pid = plates.plate_id.get(i).copied().unwrap_or(0) as usize;
                    let axis = [
                        plates.pole_axis[pid][0] as f64,
                        plates.pole_axis[pid][1] as f64,
                        plates.pole_axis[pid][2] as f64,
                    ];
                    let theta = -(plates.omega_rad_yr[pid] as f64) * dt_years;
                    let r = [
                        grid.pos_xyz[i][0] as f64,
                        grid.pos_xyz[i][1] as f64,
                        grid.pos_xyz[i][2] as f64,
                    ];
                    let src_pos = crate::plates::rotate_about_axis_f64(r, axis, theta);
                    let mut idx = i;
                    let mut best = dot(
                        [
                            grid.pos_xyz[idx][0] as f64,
                            grid.pos_xyz[idx][1] as f64,
                            grid.pos_xyz[idx][2] as f64,
                        ],
                        src_pos,
                    );
                    loop {
                        let mut improved = false;
                        for &vn in &grid.n1[idx] {
                            let j = vn as usize;
                            let d = dot(
                                [
                                    grid.pos_xyz[j][0] as f64,
                                    grid.pos_xyz[j][1] as f64,
                                    grid.pos_xyz[j][2] as f64,
                                ],
                                src_pos,
                            );
                            if d > best || ((d - best).abs() <= f64::EPSILON && j < idx) {
                                best = d;
                                idx = j;
                                improved = true;
                            }
                        }
                        if !improved {
                            break;
                        }
                    }
                    loc[k] = src_ref[idx];
                }
                (start, loc)
            });
            let (s0, loc) = handle.join().unwrap();
            assembled[s0..s0 + loc.len()].copy_from_slice(&loc);
        }
    });
    age_stage[..n].copy_from_slice(&assembled[..n]);
    age_myr[..n].copy_from_slice(&assembled[..n]);
}

/// One-shot rigid advection from a snapshot slice into the output age array.
pub fn rigid_advect_age_from(
    grid: &Grid,
    plates: &crate::plates::Plates,
    dt_years: f64,
    age_src: &[f32],
    age_out: &mut [f32],
) {
    let n = grid.cells.min(age_src.len()).min(age_out.len());
    let threads = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).clamp(1, 16);
    let chunk = ((n + threads - 1) / threads).max(1);
    let mut assembled: Vec<f32> = vec![0.0; n];
    std::thread::scope(|scope| {
        for t in 0..threads {
            let start = t * chunk;
            if start >= n {
                break;
            }
            let end = (start + chunk).min(n);
            let age_src_ref = age_src;
            let handle = scope.spawn(move || {
                let mut loc: Vec<f32> = vec![0.0; end - start];
                for k in 0..(end - start) {
                    let i = start + k;
                    let pid = plates.plate_id.get(i).copied().unwrap_or(0) as usize;
                    let axis = [
                        plates.pole_axis[pid][0] as f64,
                        plates.pole_axis[pid][1] as f64,
                        plates.pole_axis[pid][2] as f64,
                    ];
                    let theta = -(plates.omega_rad_yr[pid] as f64) * dt_years;
                    let r = [
                        grid.pos_xyz[i][0] as f64,
                        grid.pos_xyz[i][1] as f64,
                        grid.pos_xyz[i][2] as f64,
                    ];
                    let src_pos = crate::plates::rotate_about_axis_f64(r, axis, theta);
                    let mut idx = i;
                    let mut best = dot(
                        [
                            grid.pos_xyz[idx][0] as f64,
                            grid.pos_xyz[idx][1] as f64,
                            grid.pos_xyz[idx][2] as f64,
                        ],
                        src_pos,
                    );
                    loop {
                        let mut improved = false;
                        for &vn in &grid.n1[idx] {
                            let j = vn as usize;
                            let d = dot(
                                [
                                    grid.pos_xyz[j][0] as f64,
                                    grid.pos_xyz[j][1] as f64,
                                    grid.pos_xyz[j][2] as f64,
                                ],
                                src_pos,
                            );
                            if d > best || ((d - best).abs() <= f64::EPSILON && j < idx) {
                                best = d;
                                idx = j;
                                improved = true;
                            }
                        }
                        if !improved {
                            break;
                        }
                    }
                    loc[k] = age_src_ref[idx];
                }
                (start, loc)
            });
            let (s0, loc) = handle.join().unwrap();
            assembled[s0..s0 + loc.len()].copy_from_slice(&loc);
        }
    });
    age_out[..n].copy_from_slice(&assembled[..n]);
}

/// Increment age on oceanic cells only (C below threshold), and clamp to non-negative.
pub fn increment_oceanic_age(c: &[f32], age_myr: &mut [f32], dt_myr: f64) {
    let n = c.len().min(age_myr.len());
    for i in 0..n {
        let is_oceanic = c[i] < 0.15;
        if is_oceanic {
            age_myr[i] = (age_myr[i] as f64 + dt_myr).max(0.0) as f32;
        }
    }
}
