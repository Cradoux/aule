//! Continental rifting and passive margins (MVP, deterministic).

use crate::{
    boundaries::{Boundaries, EdgeClass},
    geo,
    grid::Grid,
};

const KM: f64 = 1000.0;
const RADIUS_M: f64 = 6_371_000.0;

/// Parameters controlling bilateral rift thinning and oceanization.
#[derive(Clone, Copy, Debug)]
pub struct RiftingParams {
    /// Minimum continental fraction on either side to allow rifting
    pub c_rift_min: f32,
    /// Minimum opening speed (m/yr) to be considered (noise threshold)
    pub v_open_min_m_per_yr: f32,
    /// Core half-width (km)
    pub w_core_km: f32,
    /// Margin taper width (km) beyond core
    pub w_taper_km: f32,
    /// Thinning efficiency factor (fraction of separation into crustal loss)
    pub k_thin: f32,
    /// Airy-like subsidence scaling factor
    pub alpha_subs: f32,
    /// Oceanization C threshold
    pub ocean_thresh: f32,
    /// Continental fraction decay rate toward ocean inside core (per Myr)
    pub k_c_oceanize: f32,
    /// Reset oceanic age when core becomes oceanic
    pub reset_age_on_core: bool,
    /// Enable shoulder uplift bands
    pub enable_shoulder: bool,
    /// Shoulder band width (km) beyond core (per side)
    pub w_bulge_km: f32,
    /// Shoulder uplift coefficient (fraction of thinning)
    pub beta_shoulder: f32,
    /// Couple added mass to flexure (not implemented here)
    pub couple_flexure: bool,
    /// Minimum crust thickness (m)
    pub thc_min_m: f32,
    /// Maximum crust thickness (m)
    pub thc_max_m: f32,
}

/// Stats summary for HUD/logs.
#[derive(Default, Clone, Copy, Debug)]
pub struct RiftingStats {
    /// Divergent continental edges
    pub edges_rift: usize,
    /// Cells within core bands (either side)
    pub cells_core: usize,
    /// Cells within margin taper bands (either side)
    pub cells_margin: usize,
    /// Area-weighted mean Δth_c (m), negative for thinning
    pub dthc_mean_m: f32,
    /// Min Δth_c (most negative) (m)
    pub dthc_min_m: f32,
    /// Area-weighted mean subsidence (m, positive down)
    pub subs_mean_m: f32,
    /// Max subsidence (m)
    pub subs_max_m: f32,
    /// Cells driven below oceanization threshold
    pub oceanized_cells: usize,
    /// Cells with age reset due to core oceanization
    pub age_reset_cells: usize,
}

/// Apply rifting to continental divergent zones.
#[allow(clippy::too_many_arguments)]
pub fn apply_rifting(
    grid: &Grid,
    boundaries: &Boundaries,
    plate_id: &[u16],
    _vel_m_per_yr: &[[f32; 3]],
    c: &mut [f32],
    th_c_m: &mut [f32],
    age_ocean_myr: &mut [f32],
    depth_m: &mut [f32],
    area_m2: &[f32],
    p: &RiftingParams,
    dt_myr: f64,
) -> RiftingStats {
    let n = grid.cells;
    if n == 0 {
        return RiftingStats::default();
    }

    // Identify divergent continental edges and collect seeds per side
    let mut seeds_left: Vec<u32> = Vec::new();
    let mut seeds_right: Vec<u32> = Vec::new();
    let mut edges_rift = 0usize;
    let mut vd_acc = 0.0f64;
    let mut vd_cnt = 0usize;
    for ek in &boundaries.edge_kin {
        if ek.class != EdgeClass::Divergent {
            continue;
        }
        let (u, v) = (ek.u as usize, ek.v as usize);
        let cl = c[u];
        let cr = c[v];
        if cl.max(cr) < p.c_rift_min {
            continue;
        }
        let vd = (ek.n_m_per_yr as f64).max(0.0);
        if vd < (p.v_open_min_m_per_yr as f64) {
            continue;
        }
        // Store as left/right deterministically (u left, v right) for dual distances
        seeds_left.push(ek.u);
        seeds_right.push(ek.v);
        vd_acc += vd;
        vd_cnt += 1;
        edges_rift += 1;
    }
    if edges_rift == 0 {
        return RiftingStats::default();
    }
    let vd_mean = if vd_cnt > 0 { vd_acc / (vd_cnt as f64) } else { 0.0 };

    // Precompute unit positions
    let mut rhat: Vec<[f64; 3]> = Vec::with_capacity(n);
    for r in &grid.pos_xyz {
        rhat.push(geo::normalize([r[0] as f64, r[1] as f64, r[2] as f64]));
    }

    // Dual-plate Dijkstra distances constrained by plate id of the first seed on each side
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

    if let Some(&s0) = seeds_left.first() {
        let pid = plate_id[s0 as usize];
        run_dijkstra_plate(&mut dist_left, &seeds_left, pid);
    }
    if let Some(&s0) = seeds_right.first() {
        let pid = plate_id[s0 as usize];
        run_dijkstra_plate(&mut dist_right, &seeds_right, pid);
    }

    // Weights and edits
    let w_core_m = (p.w_core_km as f64) * KM;
    let w_taper_m = (p.w_taper_km as f64) * KM;
    let w_bulge_m = (p.w_bulge_km as f64) * KM;
    let r_thin = (p.k_thin as f64) * (vd_mean * 1.0e6);
    let dt = dt_myr.max(0.0);

    let mut dthc_area = 0.0f64;
    let mut dthc_min = 0.0f64;
    let mut subs_area = 0.0f64;
    let mut subs_max = 0.0f64;
    let mut cells_core = 0usize;
    let mut cells_margin = 0usize;
    let mut oceanized_cells = 0usize;
    let mut age_reset_cells = 0usize;

    for i in 0..n {
        let dl = dist_left[i];
        let dr = dist_right[i];
        let side_w = |d: f64| -> (f64, bool, bool, bool) {
            if !d.is_finite() {
                return (0.0, false, false, false);
            }
            if d <= w_core_m {
                // cos^2 profile
                let x = std::f64::consts::PI * d / (2.0 * w_core_m);
                return ((x.cos()).powi(2), true, false, false);
            }
            if d <= w_core_m + w_taper_m {
                return (1.0 - (d - w_core_m) / w_taper_m, false, true, false);
            }
            if d <= w_core_m + w_bulge_m {
                return (1.0 - (d - w_core_m) / w_bulge_m, false, false, true);
            }
            (0.0, false, false, false)
        };
        let (wl, in_core_l, in_taper_l, in_bulge_l) = side_w(dl);
        let (wr, in_core_r, in_taper_r, in_bulge_r) = side_w(dr);
        let w_sum = (wl + wr).clamp(0.0, 1.0);
        if in_core_l || in_core_r {
            cells_core += 1;
        }
        if in_taper_l || in_taper_r {
            cells_margin += 1;
        }

        let mut dthc = 0.0f64;
        let mut subs = 0.0f64;
        if w_sum > 0.0 {
            // thinning distributed by side weights
            let w_side = w_sum; // symmetric application
            dthc -= r_thin * dt * w_side;
            let th_new = (th_c_m[i] as f64 + dthc).clamp(p.thc_min_m as f64, p.thc_max_m as f64);
            let removed = (th_c_m[i] as f64 - th_new).max(0.0);
            th_c_m[i] = th_new as f32;
            // Airy-like subsidence; positive down
            subs += (p.alpha_subs as f64) * removed;
            depth_m[i] = (depth_m[i] + subs as f32).clamp(-8000.0, 8000.0);
            dthc_area += dthc * (area_m2[i] as f64);
            subs_area += subs * (area_m2[i] as f64);
            dthc_min = dthc_min.min(dthc);
            subs_max = subs_max.max(subs);

            // Oceanization at core
            if w_sum >= 0.9 {
                let c_before = c[i] as f64;
                let c_after = (c_before - (p.k_c_oceanize as f64) * dt).max(0.0);
                c[i] = c_after as f32;
                if c_before >= (p.ocean_thresh as f64) && c_after < (p.ocean_thresh as f64) {
                    oceanized_cells += 1;
                    if p.reset_age_on_core {
                        age_ocean_myr[i] = 0.0;
                        age_reset_cells += 1;
                    }
                }
            }
        }
        // Shoulder uplift bands
        if p.enable_shoulder && (in_bulge_l || in_bulge_r) {
            let w_b: f64 = (if in_bulge_l { 1.0 } else { 0.0 }) + (if in_bulge_r { 1.0 } else { 0.0 });
            let uplift = (p.beta_shoulder as f64) * r_thin * dt * w_b.min(1.0);
            depth_m[i] = (depth_m[i] - uplift as f32).clamp(-8000.0, 8000.0);
        }
    }

    let area_tot: f64 = area_m2.iter().map(|&a| a as f64).sum();
    let stats = RiftingStats {
        edges_rift,
        cells_core,
        cells_margin,
        dthc_mean_m: if area_tot > 0.0 { (dthc_area / area_tot) as f32 } else { 0.0 },
        dthc_min_m: dthc_min as f32,
        subs_mean_m: if area_tot > 0.0 { (subs_area / area_tot) as f32 } else { 0.0 },
        subs_max_m: subs_max as f32,
        oceanized_cells,
        age_reset_cells,
    };
    println!(
        "[rifting] edges={} core={} margin={} oceanized={} age0={} | dthc mean/min={:.1}/{:.1} m subs mean/max={:.1}/{:.1} m (dt={:.1} Myr)",
        stats.edges_rift, stats.cells_core, stats.cells_margin, stats.oceanized_cells, stats.age_reset_cells, stats.dthc_mean_m, stats.dthc_min_m, stats.subs_mean_m, stats.subs_max_m, dt_myr
    );
    stats
}
