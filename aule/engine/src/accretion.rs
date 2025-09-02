//! Arc/terrane accretion at O–C convergent margins (MVP, deterministic).
//!
//! Notes:
//! - Uses boolean arc/trench masks (no continuous distances) so growth is uniform within bands.
//! - `beta_arc` applies uplift proportional to thickness change directly onto `depth_m`.
//! - Terrane docking transfers `C` locally and is highly schematic; no momentum/force coupling.

use crate::{
    boundaries::{Boundaries, EdgeClass},
    grid::Grid,
    subduction::SubductionMasks,
};

/// Parameters controlling accretion growth.
#[derive(Clone, Copy, Debug)]
pub struct AccretionParams {
    /// Thickness growth scaling at arc: r = k_arc * (Vn * 1e6) [m/Myr]
    pub k_arc: f32,
    /// Obliquity exponent for cos^gamma damping
    pub gamma_obliquity: f32,
    /// Uplift-to-thickness proportionality (applied to depth_m as negative)
    pub beta_arc: f32,
    /// Fractional C gain per reference thickness at arc
    pub alpha_arc: f32,
    /// Fractional C gain per reference thickness in forearc wedge
    pub alpha_forearc: f32,
    /// Continental C threshold
    pub c_min_continent: f32,
    /// Minimum/maximum crust thickness clamps (m)
    pub thc_min_m: f32,
    /// Maximum crust thickness (m)
    pub thc_max_m: f32,
    /// Enable terrane docking transfer
    pub enable_docking: bool,
    /// Terrane C minimum on subducting side to count as micro-continent
    pub c_terrane_min: f32,
    /// Docking distance threshold from trench (km)
    pub d_dock_km: f32,
    /// Minimum Vn for docking (m/yr)
    pub vn_min_m_per_yr: f32,
    /// Docked fraction per event
    pub tau_dock: f32,
    /// Couple added mass to flexure (not implemented here)
    pub couple_flexure: bool,
}

/// Summary stats for HUD/logs.
#[derive(Default, Clone, Copy, Debug)]
pub struct AccretionStats {
    /// Number of O–C edges considered
    pub edges_oc: usize,
    /// Cells touched in arc band
    pub cells_arc: usize,
    /// Cells touched in forearc band
    pub cells_forearc: usize,
    /// Area-weighted mean Δth_c (m)
    pub dthc_mean_m: f32,
    /// Max Δth_c (m)
    pub dthc_max_m: f32,
    /// Area-weighted mean ΔC
    pub d_c_mean: f32,
    /// Max ΔC
    pub d_c_max: f32,
    /// Docking event count
    pub dock_events: usize,
}

/// Apply accretion on overriding side at O–C convergent margins.
#[allow(clippy::too_many_arguments)]
pub fn apply_oc_accretion(
    grid: &Grid,
    subduction: &SubductionMasks,
    boundaries: &Boundaries,
    _vel_m_per_yr: &[[f32; 3]],
    c: &mut [f32],
    th_c_m: &mut [f32],
    depth_m: &mut [f32],
    area_m2: &[f32],
    p: &AccretionParams,
    dt_myr: f64,
) -> AccretionStats {
    let n = grid.cells;
    if n == 0 {
        return AccretionStats::default();
    }

    // Identify O–C edges via C thresholds per side using boundary edge list.
    let mut edges_oc = 0usize;
    let mut vn_mean: f64 = 0.0;
    let mut vn_cnt: usize = 0;
    for ek in &boundaries.edge_kin {
        if ek.class != EdgeClass::Convergent {
            continue;
        }
        let (u, v) = (ek.u as usize, ek.v as usize);
        let cl = c[u];
        let cr = c[v];
        let is_cont_l = cl >= p.c_min_continent;
        let is_cont_r = cr >= p.c_min_continent;
        if is_cont_l == is_cont_r {
            continue;
        }
        let vn = (-ek.n_m_per_yr).max(0.0) as f64;
        if vn <= 0.0 {
            continue;
        }
        vn_mean += vn;
        vn_cnt += 1;
        edges_oc += 1;
    }
    if edges_oc == 0 {
        return AccretionStats::default();
    }
    let vn_ref = if vn_cnt > 0 { vn_mean / (vn_cnt as f64) } else { 0.0 };

    // Spatial weights from masks: arc/backarc on overriding, trench–arc forearc band.
    // We approximate forearc by cells that are not arc but within trench OR arc masks proximity.
    let mut dthc_area = 0.0f64;
    let mut dthc_max = 0.0f64;
    let mut d_c_area = 0.0f64;
    let mut d_c_max = 0.0f64;
    let mut cells_arc = 0usize;
    let mut cells_forearc = 0usize;
    let dt = dt_myr.max(0.0);
    const H_REF: f64 = 35_000.0;

    // Gaussian width for arc from subduction params: sigma = half_width/2. Use boolean mask as proxy.
    // Without explicit distances to arc center here, use a uniform weight inside the arc mask.
    let r_base = (p.k_arc as f64) * (vn_ref * 1.0e6);

    for i in 0..n {
        let mut dthc = 0.0f64;
        let mut dcf = 0.0f64;
        if subduction.arc[i] {
            cells_arc += 1;
            let w_arc = 1.0f64; // proxy weight
            dthc += (p.beta_arc as f64) * r_base * dt * w_arc;
            dcf += (p.alpha_arc as f64) * (r_base / H_REF) * dt * w_arc;
        }
        // Forearc: between trench and arc; approximate as trench OR arc neighbor cells not in arc
        if subduction.trench[i] && !subduction.arc[i] {
            cells_forearc += 1;
            let w_fore = 1.0f64; // proxy weight
            dcf += (p.alpha_forearc as f64) * ((vn_ref * 1.0e6) / H_REF) * dt * w_fore;
        }
        if dthc != 0.0 || dcf != 0.0 {
            // Per-operator safety: clamp per-step thickness and uplift; flag units issues
            let cap_th = 200.0f64;
            let cap_uplift = 200.0f64;
            if dthc.abs() > 1000.0 {
                println!("[UNITS_BUG] accretion Δth_c raw {:.1} m (>1000 m)", dthc);
            }
            let dthc_c = dthc.clamp(-cap_th, cap_th);
            let uplift = dthc_c; // beta already applied into dthc
            let th_new = (th_c_m[i] as f64 + dthc_c).clamp(p.thc_min_m as f64, p.thc_max_m as f64);
            let dc_new = (c[i] as f64 + dcf).clamp(0.0, 1.0);
            th_c_m[i] = th_new as f32;
            c[i] = dc_new as f32;
            let uplift_c = uplift.clamp(-cap_uplift, cap_uplift);
            depth_m[i] = (depth_m[i] - uplift_c as f32).clamp(-8000.0, 8000.0);
            dthc_area += dthc * (area_m2[i] as f64);
            d_c_area += dcf * (area_m2[i] as f64);
            dthc_max = dthc_max.max(dthc);
            d_c_max = d_c_max.max(dcf);
        }
    }

    // Docking MVP: count events near trench where subducting side has high C
    let mut dock_events = 0usize;
    if p.enable_docking {
        for ek in &boundaries.edge_kin {
            if ek.class != EdgeClass::Convergent {
                continue;
            }
            let (u, v) = (ek.u as usize, ek.v as usize);
            let cl = c[u];
            let cr = c[v];
            let is_cont_l = cl >= p.c_terrane_min;
            let is_cont_r = cr >= p.c_terrane_min;
            let vn = (-ek.n_m_per_yr).max(0.0);
            if vn < p.vn_min_m_per_yr {
                continue;
            }
            // One side must be terrane and near trench mask; approximate using the cells themselves.
            if subduction.trench[u] && is_cont_l && !subduction.trench[v] {
                let tau = p.tau_dock.clamp(0.0, 0.05);
                let give = (c[u] * tau).min(1.0 - c[v]);
                c[u] -= give;
                c[v] += give;
                dock_events += 1;
            } else if subduction.trench[v] && is_cont_r && !subduction.trench[u] {
                let tau = p.tau_dock.clamp(0.0, 0.05);
                let give = (c[v] * tau).min(1.0 - c[u]);
                c[v] -= give;
                c[u] += give;
                dock_events += 1;
            }
        }
    }

    let area_tot: f64 = area_m2.iter().map(|&a| a as f64).sum();
    let stats = AccretionStats {
        edges_oc,
        cells_arc,
        cells_forearc,
        dthc_mean_m: if area_tot > 0.0 { (dthc_area / area_tot) as f32 } else { 0.0 },
        dthc_max_m: dthc_max as f32,
        d_c_mean: if area_tot > 0.0 { (d_c_area / area_tot) as f32 } else { 0.0 },
        d_c_max: d_c_max as f32,
        dock_events,
    };
    println!(
        "[accretion] OC edges={} arc={} forearc={} | dthc mean/max={:.1}/{:.1} m dC mean/max={:.4}/{:.4} | dock={} (dt={:.1} Myr)",
        stats.edges_oc, stats.cells_arc, stats.cells_forearc, stats.dthc_mean_m, stats.dthc_max_m, stats.d_c_mean, stats.d_c_max, stats.dock_events, dt_myr
    );
    stats
}
