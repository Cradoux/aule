//! Plate boundary classification (T-040).
//! Classifies edges between different plates as divergent, convergent, or transform.

use crate::grid::Grid;

/// Summary counts of boundary classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryStats {
    /// Divergent boundary edges count.
    pub divergent: u32,
    /// Convergent boundary edges count.
    pub convergent: u32,
    /// Transform boundary edges count.
    pub transform: u32,
}

impl BoundaryStats {
    fn zero() -> Self {
        Self { divergent: 0, convergent: 0, transform: 0 }
    }
}

/// Boundary edge class.
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClass {
    /// Divergent boundary
    Divergent = 1,
    /// Convergent boundary
    Convergent = 2,
    /// Transform boundary
    Transform = 3,
}

impl From<u8> for EdgeClass {
    fn from(v: u8) -> Self {
        match v {
            1 => Self::Divergent,
            2 => Self::Convergent,
            3 => Self::Transform,
            _ => Self::Transform,
        }
    }
}

/// Per-edge stored kinematics.
pub struct EdgeKin {
    /// First cell index (u < v)
    pub u: u32,
    /// Second cell index
    pub v: u32,
    /// Edge class
    pub class: EdgeClass,
    /// Across-boundary unit normal in world coords
    pub n_hat: [f32; 3],
    /// Along-boundary unit tangent in world coords
    pub t_hat: [f32; 3],
    /// Normal component of relative velocity (m/yr), signed (positive opening)
    pub n_m_per_yr: f32,
    /// Tangential component magnitude of relative velocity (m/yr)
    pub t_m_per_yr: f32,
}

/// Classified boundaries.
pub struct Boundaries {
    /// Per-cell bitmask: bit0=div, bit1=conv, bit2=trans.
    pub b: Vec<u8>,
    /// Optional compact undirected edge list (u < v), class in {1:div,2:conv,3:trans}.
    pub edges: Vec<(u32, u32, u8)>,
    /// Persisted per-edge kinematics
    pub edge_kin: Vec<EdgeKin>,
    /// Summary statistics.
    pub stats: BoundaryStats,
}

impl Boundaries {
    /// Classify boundaries from `plate_id` and per-cell local velocities `v_en` (east,north) in m/yr.
    /// `tau_open_m_per_yr` is the opening/closing threshold (e.g., 0.005 = 0.5 cm/yr).
    pub fn classify(
        grid: &Grid,
        plate_id: &[u16],
        v_en: &[[f32; 2]],
        tau_open_m_per_yr: f64,
    ) -> Self {
        assert_eq!(plate_id.len(), grid.cells);
        assert_eq!(v_en.len(), grid.cells);

        let tau = tau_open_m_per_yr.max(0.0);

        // Precompute world positions and local EN bases per cell.
        let mut pos: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
        let mut east: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
        let mut north: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
        for i in 0..grid.cells {
            let r =
                [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
            pos.push(r);
            let e = normalize(cross([0.0, 0.0, 1.0], r));
            let n = normalize(cross(r, e));
            east.push(e);
            north.push(n);
        }

        let mut b = vec![0u8; grid.cells];
        let mut edges: Vec<(u32, u32, u8)> = Vec::new();
        let mut stats = BoundaryStats::zero();
        let mut edge_kin: Vec<EdgeKin> = Vec::new();

        for u in 0..grid.cells as u32 {
            for &v in &grid.n1[u as usize] {
                if v <= u {
                    continue;
                }
                if plate_id[u as usize] == plate_id[v as usize] {
                    continue;
                }

                let ru = pos[u as usize];
                let rv = pos[v as usize];
                // Midpoint great-circle geometry
                let rm = normalize([ru[0] + rv[0], ru[1] + rv[1], ru[2] + rv[2]]);
                let g = normalize(cross(ru, rv));
                let t_hat = normalize(cross(g, rm)); // along-boundary
                let mut n_hat = normalize(cross(rm, t_hat)); // across-boundary
                if dot(n_hat, [rv[0] - ru[0], rv[1] - ru[1], rv[2] - ru[2]]) < 0.0 {
                    n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                }

                // Convert local EN velocities to world vectors
                let vu = [
                    (v_en[u as usize][0] as f64) * east[u as usize][0]
                        + (v_en[u as usize][1] as f64) * north[u as usize][0],
                    (v_en[u as usize][0] as f64) * east[u as usize][1]
                        + (v_en[u as usize][1] as f64) * north[u as usize][1],
                    (v_en[u as usize][0] as f64) * east[u as usize][2]
                        + (v_en[u as usize][1] as f64) * north[u as usize][2],
                ];
                let vv = [
                    (v_en[v as usize][0] as f64) * east[v as usize][0]
                        + (v_en[v as usize][1] as f64) * north[v as usize][0],
                    (v_en[v as usize][0] as f64) * east[v as usize][1]
                        + (v_en[v as usize][1] as f64) * north[v as usize][1],
                    (v_en[v as usize][0] as f64) * east[v as usize][2]
                        + (v_en[v as usize][1] as f64) * north[v as usize][2],
                ];
                let dv = [vu[0] - vv[0], vu[1] - vv[1], vu[2] - vv[2]];
                let n = dot(dv, n_hat);
                let t = dot(dv, t_hat).abs();

                let class = if n > tau {
                    1u8 // divergent
                } else if n < -tau {
                    2u8 // convergent
                } else if t > n.abs() {
                    3u8 // transform
                } else {
                    0u8 // below threshold; ignore
                };

                if class != 0 {
                    // per-cell bits: bit0=div, bit1=conv, bit2=trans
                    match class {
                        1 => {
                            b[u as usize] |= 1;
                            b[v as usize] |= 1;
                            stats.divergent += 1;
                        }
                        2 => {
                            b[u as usize] |= 1 << 1;
                            b[v as usize] |= 1 << 1;
                            stats.convergent += 1;
                        }
                        3 => {
                            b[u as usize] |= 1 << 2;
                            b[v as usize] |= 1 << 2;
                            stats.transform += 1;
                        }
                        _ => {}
                    }
                    edges.push((u, v, class));
                    edge_kin.push(EdgeKin {
                        u,
                        v,
                        class: EdgeClass::from(class),
                        n_hat: [n_hat[0] as f32, n_hat[1] as f32, n_hat[2] as f32],
                        t_hat: [t_hat[0] as f32, t_hat[1] as f32, t_hat[2] as f32],
                        n_m_per_yr: n as f32,
                        t_m_per_yr: t as f32,
                    });
                }
            }
        }

        Self { b, edges, edge_kin, stats }
    }
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
