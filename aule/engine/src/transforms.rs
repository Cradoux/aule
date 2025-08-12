//! Transform pull-apart vs restraining bands and bathymetry edits (CPU, deterministic).

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::{boundaries::Boundaries, geo, grid::Grid};

const KM: f64 = 1000.0;
const RADIUS_M: f64 = 6_371_000.0;

/// Boolean per-cell masks for transform features.
#[derive(Default)]
pub struct TransformMasks {
    /// Cells in pull-apart (extensional) basins.
    pub pull_apart: Vec<bool>,
    /// Cells in restraining (compressional) uplifts.
    pub restraining: Vec<bool>,
}

/// Parameters controlling transform detection and band strengths.
#[derive(Clone, Copy)]
pub struct TransformParams {
    /// Opening/closing threshold (m/yr) to consider near-pure transform.
    pub tau_open_m_per_yr: f64,
    /// Minimum tangential speed (m/yr) to deem transform active.
    pub min_tangential_m_per_yr: f64,
    /// Half-width of transform band (km) measured from the fault trace.
    pub basin_half_width_km: f64,
    /// Uplift (negative, shallower) applied in restraining zones (m).
    pub ridge_like_uplift_m: f32,
    /// Deepening (positive, deeper) applied in pull-apart basins (m).
    pub basin_deepen_m: f32,
}

/// Summary statistics for transform band sizes.
pub struct TransformStats {
    /// Number of pull-apart cells.
    pub pull_apart_cells: u32,
    /// Number of restraining cells.
    pub restraining_cells: u32,
}

/// Detect transform segments and apply simple bathymetry edits.
pub fn apply_transforms(
    grid: &Grid,
    boundaries: &Boundaries,
    v_en: &[[f32; 2]],
    depth_m: &mut [f32],
    params: TransformParams,
) -> (TransformMasks, TransformStats) {
    assert_eq!(v_en.len(), grid.cells);
    assert_eq!(depth_m.len(), grid.cells);

    // Precompute unit positions and world velocities
    let mut rhat: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
    for p in &grid.pos_xyz {
        rhat.push(geo::normalize([p[0] as f64, p[1] as f64, p[2] as f64]));
    }
    let vel_w: Vec<[f64; 3]> =
        rhat.iter().enumerate().map(|(i, &rh)| geo::en_to_world(rh, v_en[i])).collect();

    // Seeds for pull-apart vs restraining bands
    let mut seeds_pull: Vec<u32> = Vec::new();
    let mut seeds_rest: Vec<u32> = Vec::new();

    for &(u, v, class) in &boundaries.edges {
        if class != 3 {
            continue; // transform only
        }
        let u_idx = u as usize;
        let v_idx = v as usize;
        let rh_u = rhat[u_idx];
        let rh_v = rhat[v_idx];
        // Tangent from u to v at midpoint direction
        let t_dir = geo::normalize(geo::cross(rh_u, rh_v));
        // Relative velocity in world coords
        let dv = [
            vel_w[v_idx][0] - vel_w[u_idx][0],
            vel_w[v_idx][1] - vel_w[u_idx][1],
            vel_w[v_idx][2] - vel_w[u_idx][2],
        ];
        let n_hat = geo::normalize(geo::cross(
            t_dir,
            geo::normalize([
                (rh_u[0] + rh_v[0]) * 0.5,
                (rh_u[1] + rh_v[1]) * 0.5,
                (rh_u[2] + rh_v[2]) * 0.5,
            ]),
        ));
        let t_comp = geo::dot(dv, t_dir).abs();
        let n_comp = geo::dot(dv, n_hat).abs();
        if t_comp >= params.min_tangential_m_per_yr && n_comp <= params.tau_open_m_per_yr {
            // Shear sense using signed projection onto t_dir from u->v
            let shear_sign = geo::dot(dv, t_dir);
            if shear_sign >= 0.0 {
                seeds_pull.push(u);
                seeds_pull.push(v);
            } else {
                seeds_rest.push(u);
                seeds_rest.push(v);
            }
        }
    }

    #[derive(Copy, Clone)]
    struct QItem {
        d: f64,
        c: u32,
    }
    impl Eq for QItem {}
    impl PartialEq for QItem {
        fn eq(&self, o: &Self) -> bool {
            self.d.eq(&o.d) && self.c == o.c
        }
    }
    impl PartialOrd for QItem {
        fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
            Some(self.cmp(o))
        }
    }
    impl Ord for QItem {
        fn cmp(&self, o: &Self) -> Ordering {
            match self.d.partial_cmp(&o.d) {
                Some(Ordering::Less) => Ordering::Greater,
                Some(Ordering::Greater) => Ordering::Less,
                Some(Ordering::Equal) => self.c.cmp(&o.c).reverse(),
                None => Ordering::Equal,
            }
        }
    }

    // Distance fields
    let mut dist_pull: Vec<f64> = vec![f64::INFINITY; grid.cells];
    let mut dist_rest: Vec<f64> = vec![f64::INFINITY; grid.cells];

    let run_dijkstra = |dist: &mut [f64], seeds: &[u32]| {
        let mut h: BinaryHeap<QItem> = BinaryHeap::new();
        for &s in seeds {
            dist[s as usize] = 0.0;
            h.push(QItem { d: 0.0, c: s });
        }
        while let Some(QItem { d, c }) = h.pop() {
            let u = c as usize;
            if d > dist[u] {
                continue;
            }
            for &vn in &grid.n1[u] {
                let v = vn as usize;
                let s_m = geo::great_circle_arc_len_m(rhat[u], rhat[v], RADIUS_M);
                let nd = d + s_m;
                if nd < dist[v] {
                    dist[v] = nd;
                    h.push(QItem { d: nd, c: v as u32 });
                }
            }
        }
    };

    if !seeds_pull.is_empty() {
        run_dijkstra(&mut dist_pull, &seeds_pull);
    }
    if !seeds_rest.is_empty() {
        run_dijkstra(&mut dist_rest, &seeds_rest);
    }

    // Build masks
    let hw_m = params.basin_half_width_km * KM;
    let mut masks = TransformMasks {
        pull_apart: vec![false; grid.cells],
        restraining: vec![false; grid.cells],
    };
    let mut stats = TransformStats { pull_apart_cells: 0, restraining_cells: 0 };
    for i in 0..grid.cells {
        if dist_pull[i].is_finite() && dist_pull[i] <= hw_m {
            masks.pull_apart[i] = true;
            stats.pull_apart_cells += 1;
        }
        if dist_rest[i].is_finite() && dist_rest[i] <= hw_m {
            masks.restraining[i] = true;
            stats.restraining_cells += 1;
        }
    }

    // Apply bathy edits where not already changed by caller (viewer can enforce priority externally)
    for (i, d) in depth_m.iter_mut().enumerate() {
        if masks.pull_apart[i] {
            *d += params.basin_deepen_m;
        } else if masks.restraining[i] {
            *d += params.ridge_like_uplift_m;
        }
    }

    (masks, stats)
}
