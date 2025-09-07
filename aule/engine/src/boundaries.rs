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

        // Use precomputed world positions and local EN bases (from Grid precompute)
        let mut pos: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
        for i in 0..grid.cells {
            pos.push([
                grid.pos_xyz[i][0] as f64,
                grid.pos_xyz[i][1] as f64,
                grid.pos_xyz[i][2] as f64,
            ]);
        }

        // Parallel classification by chunking u-range; accumulate edges then set b from edges
        let mut edges_all: Vec<(u32, u32, u8)> = Vec::new();
        let mut kin_all: Vec<EdgeKin> = Vec::new();
        let mut stats = BoundaryStats::zero();
        {
            use std::thread;
            let n_cells = grid.cells as u32;
            let threads =
                thread::available_parallelism().map(|n| n.get()).unwrap_or(1).clamp(1, 16);
            let chunk = ((n_cells as usize + threads - 1) / threads).max(1);
            let mut handles = Vec::new();
            for t in 0..threads {
                let start = (t * chunk) as u32;
                if start >= n_cells {
                    break;
                }
                let end = ((t + 1) * chunk) as u32;
                // Capture shared references
                let pos_t = pos.clone();
                let n1_t = grid.n1.clone();
                let east_t = grid.east_hat.clone();
                let north_t = grid.north_hat.clone();
                let v_en_t = v_en.to_vec();
                let pid_t = plate_id.to_vec();
                handles.push(thread::spawn(move || {
                    let mut loc_edges: Vec<(u32, u32, u8)> = Vec::new();
                    let mut loc_kin: Vec<EdgeKin> = Vec::new();
                    let mut st = BoundaryStats::zero();
                    for u in start..end.min(n_cells) {
                        for &v in &n1_t[u as usize] {
                            if v <= u {
                                continue;
                            }
                            if pid_t[u as usize] == pid_t[v as usize] {
                                continue;
                            }
                            let ru = pos_t[u as usize];
                            let rv = pos_t[v as usize];
                            let rm = crate::geo::normalize([
                                ru[0] + rv[0],
                                ru[1] + rv[1],
                                ru[2] + rv[2],
                            ]);
                            let g = crate::geo::normalize(crate::geo::cross(ru, rv));
                            let t_hat = crate::geo::normalize(crate::geo::cross(g, rm));
                            let mut n_hat = crate::geo::normalize(crate::geo::cross(rm, t_hat));
                            if crate::geo::dot(n_hat, [rv[0] - ru[0], rv[1] - ru[1], rv[2] - ru[2]])
                                < 0.0
                            {
                                n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                            }
                            let eu = [
                                east_t[u as usize][0] as f64,
                                east_t[u as usize][1] as f64,
                                east_t[u as usize][2] as f64,
                            ];
                            let nu = [
                                north_t[u as usize][0] as f64,
                                north_t[u as usize][1] as f64,
                                north_t[u as usize][2] as f64,
                            ];
                            let ev = [
                                east_t[v as usize][0] as f64,
                                east_t[v as usize][1] as f64,
                                east_t[v as usize][2] as f64,
                            ];
                            let nv = [
                                north_t[v as usize][0] as f64,
                                north_t[v as usize][1] as f64,
                                north_t[v as usize][2] as f64,
                            ];
                            let vu = [
                                (v_en_t[u as usize][0] as f64) * eu[0]
                                    + (v_en_t[u as usize][1] as f64) * nu[0],
                                (v_en_t[u as usize][0] as f64) * eu[1]
                                    + (v_en_t[u as usize][1] as f64) * nu[1],
                                (v_en_t[u as usize][0] as f64) * eu[2]
                                    + (v_en_t[u as usize][1] as f64) * nu[2],
                            ];
                            let vv = [
                                (v_en_t[v as usize][0] as f64) * ev[0]
                                    + (v_en_t[v as usize][1] as f64) * nv[0],
                                (v_en_t[v as usize][0] as f64) * ev[1]
                                    + (v_en_t[v as usize][1] as f64) * nv[1],
                                (v_en_t[v as usize][0] as f64) * ev[2]
                                    + (v_en_t[v as usize][1] as f64) * nv[2],
                            ];
                            let dv = [vu[0] - vv[0], vu[1] - vv[1], vu[2] - vv[2]];
                            let n = crate::geo::dot(dv, n_hat);
                            let t_signed = crate::geo::dot(dv, t_hat);
                            let t_abs = t_signed.abs();
                            // Classify boundary type based on velocity gradients
                            let class = if n > tau {
                                1u8 // Divergent (opening/rifting)
                            } else if n < -tau {
                                2u8 // Convergent (closing/subduction)
                            } else if t_abs > n.abs() {
                                3u8 // Transform (shearing)
                            } else {
                                0u8 // No significant activity
                            };
                            
                            if class != 0 {
                                match class {
                                    1 => {
                                        st.divergent += 1;
                                    }
                                    2 => {
                                        st.convergent += 1;
                                    }
                                    3 => {
                                        st.transform += 1;
                                    }
                                    _ => {}
                                }
                                loc_edges.push((u, v, class));
                                loc_kin.push(EdgeKin {
                                    u,
                                    v,
                                    class: EdgeClass::from(class),
                                    n_hat: [n_hat[0] as f32, n_hat[1] as f32, n_hat[2] as f32],
                                    t_hat: [t_hat[0] as f32, t_hat[1] as f32, t_hat[2] as f32],
                                    n_m_per_yr: n as f32,
                                    t_m_per_yr: t_signed as f32,
                                });
                            }
                        }
                    }
                    (loc_edges, loc_kin, st)
                }));
            }
            for h in handles {
                if let Ok((mut e, mut k, st_local)) = h.join() {
                    stats.divergent += st_local.divergent;
                    stats.convergent += st_local.convergent;
                    stats.transform += st_local.transform;
                    edges_all.append(&mut e);
                    kin_all.append(&mut k);
                }
            }
        }
        // Build per-cell bitmask from edges
        let mut b = vec![0u8; grid.cells];
        for &(u, v, class) in &edges_all {
            match class {
                1 => {
                    b[u as usize] |= 1;
                    b[v as usize] |= 1;
                }
                2 => {
                    b[u as usize] |= 1 << 1;
                    b[v as usize] |= 1 << 1;
                }
                3 => {
                    b[u as usize] |= 1 << 2;
                    b[v as usize] |= 1 << 2;
                }
                _ => {}
            }
        }
        // Deterministic ordering
        edges_all.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        kin_all.sort_unstable_by(|a, b| a.u.cmp(&b.u).then(a.v.cmp(&b.v)));

        Self { b, edges: edges_all, edge_kin: kin_all, stats }
    }
}

// Local math helpers are removed in favor of crate::geo to keep a single source of truth
