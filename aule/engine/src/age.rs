//! Steady-state oceanic age via Dijkstra from ridge seeds and age→depth mapping.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::boundaries::Boundaries;
use crate::geo;
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
