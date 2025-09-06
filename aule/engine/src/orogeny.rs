//! Collision orogeny (C–C sutures) MVP.
//!
//! Modeling notes:
//! - Thickening rate derives from normal convergence magnitude with an obliquity damping factor.
//! - Uplift is proportional to added thickness and applied as a negative `depth` delta (shallower).
//! - No explicit mass conservation between eroded material and thickening; surface processes run
//!   separately and may not balance orogenic addition. Consider coupling or running at cadence.

use crate::plates::PlateKind;
use crate::{
    boundaries::{Boundaries, EdgeClass},
    geo,
    grid::Grid,
};

const KM: f64 = 1000.0;
const RADIUS_M: f64 = 6_371_000.0;

/// Parameters controlling bilateral orogenic thickening.
#[derive(Clone, Copy, Debug)]
pub struct OrogenyParams {
    /// Minimum continental fraction on both sides to deem C–C (0..1)
    pub c_min: f32,
    /// Core half-width on each side (km)
    pub w_core_km: f32,
    /// Additional foreland taper width (km)
    pub w_taper_km: f32,
    /// Thickening efficiency (fraction of incoming normal convergence per Myr -> m/Myr)
    pub k_thick: f32,
    /// Proportionality factor from thickening to uplift (negative to shallower)
    pub beta_uplift: f32,
    /// Obliquity exponent for cos^gamma scaling
    pub gamma_obliquity: f32,
    /// If true, couple added load to flexure (not implemented in MVP)
    pub couple_flexure: bool,
}

/// Summary statistics for diagnostics and HUD readouts.
#[derive(Default, Clone, Copy, Debug)]
pub struct OrogenyStats {
    /// Number of convergent C–C edges detected
    pub edges_cc: usize,
    /// Cells influenced inside core band (either side)
    pub cells_core: usize,
    /// Cells influenced inside taper band (either side)
    pub cells_taper: usize,
    /// Area-weighted mean Δth_c (m)
    pub dthc_mean_m: f32,
    /// Max Δth_c (m)
    pub dthc_max_m: f32,
    /// Area-weighted mean uplift (m)
    pub uplift_mean_m: f32,
    /// Max uplift (m)
    pub uplift_max_m: f32,
}

/// Apply continental collision thickening and uplift.
#[allow(clippy::too_many_arguments)]
pub fn apply_cc_orogeny(
    grid: &Grid,
    boundaries: &Boundaries,
    plate_id: &[u16],
    _vel_m_per_yr: &[[f32; 3]],
    c: &[f32],
    area_m2: &[f32],
    th_c_m: &mut [f32],
    depth_m: &mut [f32],
    p: &OrogenyParams,
    dt_myr: f64,
    plates_kind: &[PlateKind],
) -> OrogenyStats {
    let n = grid.cells;
    if n == 0 {
        return OrogenyStats::default();
    }
    // Seeds per side
    let mut seeds_left: Vec<u32> = Vec::new();
    let mut seeds_right: Vec<u32> = Vec::new();
    let mut r_thick_acc: f64 = 0.0;
    let mut r_thick_cnt: usize = 0;
    let mut edges_cc = 0usize;

    for ek in &boundaries.edge_kin {
        if ek.class != EdgeClass::Convergent {
            continue;
        }
        let (ul, vr) = (ek.u as usize, ek.v as usize);
        let c_l = c.get(ul).copied().unwrap_or(0.0);
        let c_r = c.get(vr).copied().unwrap_or(0.0);
        if c_l < p.c_min || c_r < p.c_min {
            continue;
        }
        // Require both sides to be Continental by plate kind
        let pk_l = plates_kind.get(plate_id[ul] as usize).copied().unwrap_or(PlateKind::Oceanic);
        let pk_r = plates_kind.get(plate_id[vr] as usize).copied().unwrap_or(PlateKind::Oceanic);
        if !(matches!(pk_l, PlateKind::Continental) && matches!(pk_r, PlateKind::Continental)) {
            continue;
        }
        // Normal convergence magnitude (positive converge)
        let vn = (-ek.n_m_per_yr).max(0.0) as f64;
        if vn <= 0.0 {
            continue;
        }
        // Obliquity factor and width-based rate per RL-3: ḣ_c = α * v_n * sin(θ)^β / W_c
        let phi = (ek.t_m_per_yr.abs() as f64).atan2(vn);
        let sin_theta = phi.sin().abs();
        let f_obl = sin_theta.powf(p.gamma_obliquity as f64);
        let w_c_m = (p.w_core_km as f64).clamp(150.0, 300.0) * 1000.0;
        let alpha = p.k_thick as f64; // reuse k_thick as α
        let r_thick = alpha * (vn * 1.0e6) * f_obl / w_c_m; // m/Myr
        r_thick_acc += r_thick;
        r_thick_cnt += 1;
        edges_cc += 1;
        // Treat ek.u as left seed, ek.v as right seed for dual distances
        seeds_left.push(ek.u);
        seeds_right.push(ek.v);
    }

    if edges_cc == 0 {
        return OrogenyStats::default();
    }

    let r_mean = if r_thick_cnt > 0 { r_thick_acc / (r_thick_cnt as f64) } else { 0.0 };
    let w_core_m = (p.w_core_km as f64) * KM;
    let w_taper_m = (p.w_taper_km as f64) * KM;

    // Precompute rhat per cell
    let mut rhat: Vec<[f64; 3]> = Vec::with_capacity(n);
    for r in &grid.pos_xyz {
        rhat.push([r[0] as f64, r[1] as f64, r[2] as f64]);
    }

    // Dijkstra dist on each side, constrained to plate id regions using seeds' plate ids
    let mut dist_left: Vec<f64> = vec![f64::INFINITY; n];
    let mut dist_right: Vec<f64> = vec![f64::INFINITY; n];

    let run_dijkstra_plate = |dist: &mut [f64], seeds: &[u32], pid_allow: u16| {
        use std::cmp::Ordering;
        use std::collections::BinaryHeap;
        #[derive(Copy, Clone)]
        struct QI {
            d: f64,
            c: u32,
        }
        impl Eq for QI {}
        impl PartialEq for QI {
            fn eq(&self, o: &Self) -> bool {
                self.d.eq(&o.d) && self.c == o.c
            }
        }
        impl PartialOrd for QI {
            fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
                Some(self.cmp(o))
            }
        }
        impl Ord for QI {
            fn cmp(&self, o: &Self) -> Ordering {
                match self.d.partial_cmp(&o.d) {
                    Some(Ordering::Less) => Ordering::Greater,
                    Some(Ordering::Greater) => Ordering::Less,
                    Some(Ordering::Equal) => self.c.cmp(&o.c).reverse(),
                    None => Ordering::Equal,
                }
            }
        }
        let mut h: BinaryHeap<QI> = BinaryHeap::new();
        for &s in seeds {
            dist[s as usize] = 0.0;
            h.push(QI { d: 0.0, c: s });
        }
        while let Some(QI { d, c }) = h.pop() {
            let u = c as usize;
            if d > dist[u] {
                continue;
            }
            for &vn in &grid.n1[u] {
                let v = vn as usize;
                if plate_id[v] != pid_allow {
                    continue;
                }
                let s_m = geo::great_circle_arc_len_m(rhat[u], rhat[v], RADIUS_M);
                let nd = d + s_m;
                if nd < dist[v] {
                    dist[v] = nd;
                    h.push(QI { d: nd, c: v as u32 });
                }
            }
        }
    };

    if let Some(&seed0) = seeds_left.first() {
        let pid = plate_id[seed0 as usize];
        run_dijkstra_plate(&mut dist_left, &seeds_left, pid);
    }
    if let Some(&seed0) = seeds_right.first() {
        let pid = plate_id[seed0 as usize];
        run_dijkstra_plate(&mut dist_right, &seeds_right, pid);
    }

    // Apply bilateral weights
    let mut dthc_mean = 0.0f64;
    let mut dthc_max = 0.0f64;
    let mut upl_mean = 0.0f64;
    let mut upl_max = 0.0f64;
    let mut cells_core = 0usize;
    let mut cells_taper = 0usize;
    let dt = dt_myr.max(0.0);
    for i in 0..n {
        let dl = dist_left[i];
        let dr = dist_right[i];
        // compute bilateral weight sum (0..1)
        let core_m = w_core_m;
        let taper_m = w_taper_m;
        let side_w = |d: f64| -> (f64, bool, bool) {
            if d.is_finite() {
                if d <= core_m {
                    (0.5 * (1.0 + (std::f64::consts::PI * d / core_m).cos()), true, false)
                } else if d <= core_m + taper_m {
                    (1.0 - ((d - core_m) / taper_m), false, true)
                } else {
                    (0.0, false, false)
                }
            } else {
                (0.0, false, false)
            }
        };
        let (wl, in_core_l, in_taper_l) = side_w(dl);
        let (wr, in_core_r, in_taper_r) = side_w(dr);
        if in_core_l || in_core_r {
            cells_core += 1;
        }
        if in_taper_l || in_taper_r {
            cells_taper += 1;
        }
        let w_sum: f64 = (wl + wr).clamp(0.0, 1.0);
        if w_sum <= 0.0 {
            continue;
        }
        let dthc_raw = (r_mean * dt) * w_sum;
        let uplift_raw = (p.beta_uplift as f64) * dthc_raw;
        // RL-2: apply smooth saturator on rates before integrating (convert back from per-step Δ)
        let rate_thick = dthc_raw / dt.max(1e-12);
        let rate_uplift = uplift_raw / dt.max(1e-12);
        let cap_m_per_myr = 300.0f64 / dt.max(1e-12);
        let rate_thick_soft = crate::util::soft_cap_f64(rate_thick, cap_m_per_myr);
        let rate_uplift_soft = crate::util::soft_cap_f64(rate_uplift, cap_m_per_myr);
        let dthc = rate_thick_soft * dt;
        let uplift = rate_uplift_soft * dt;
        th_c_m[i] = (th_c_m[i] + dthc as f32).clamp(0.0, 70_000.0);
        depth_m[i] = (depth_m[i] - uplift as f32).clamp(-8000.0, 8000.0);
        dthc_mean += dthc * (area_m2[i] as f64);
        upl_mean += uplift * (area_m2[i] as f64);
        dthc_max = dthc_max.max(dthc);
        upl_max = upl_max.max(uplift);
    }
    let area_tot: f64 = area_m2.iter().map(|&a| a as f64).sum();
    let dthc_mean_m = if area_tot > 0.0 { (dthc_mean / area_tot) as f32 } else { 0.0 };
    let uplift_mean_m = if area_tot > 0.0 { (upl_mean / area_tot) as f32 } else { 0.0 };
    let stats = OrogenyStats {
        edges_cc,
        cells_core,
        cells_taper,
        dthc_mean_m,
        dthc_max_m: dthc_max as f32,
        uplift_mean_m,
        uplift_max_m: upl_max as f32,
    };
    println!(
        "[orogeny] CC edges={} core={} taper={} | dthc mean/max={:.1}/{:.1} m | uplift mean/max={:.1}/{:.1} m (dt={:.1} Myr)",
        stats.edges_cc, stats.cells_core, stats.cells_taper, stats.dthc_mean_m, stats.dthc_max_m, stats.uplift_mean_m, stats.uplift_max_m, dt_myr
    );
    stats
}
