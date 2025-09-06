//! Surface processes orchestrator (T-450): fluvial erosion, hillslope diffusion,
//! and sediment transport/deposition acting on `depth_m` and `sediment_m`.
//!
//! Conventions:
//! - depth_m: positive values indicate water depth below sea level; negative values indicate
//!   land elevation above sea level with magnitude in meters. Elevation is defined as
//!   `elev_m = -depth_m`.
//! - Erosion removes material and lowers elevation → increases `depth_m` (Δdepth > 0).
//! - Deposition adds material and raises elevation → decreases `depth_m` (Δdepth < 0).

use crate::grid::Grid;

/// Parameters controlling surface processes per step.
///
/// Units: distances in meters, areas in m², volumes in m³, time in years (yr).
#[derive(Clone, Copy, Debug)]
pub struct SurfaceParams {
    /// Stream-power erodibility coefficient (m^(1−m) yr^(−1)).
    pub k_stream: f32,
    /// Stream-power drainage-area exponent m (0.3…0.6).
    pub m_exp: f32,
    /// Stream-power slope exponent n (0.8…1.5).
    pub n_exp: f32,
    /// Hillslope diffusivity κ (m²/yr).
    pub k_diff: f32,
    /// Transport-capacity scale coefficient (dimensionless scale for MVP).
    pub k_tr: f32,
    /// Transport capacity drainage-area exponent p (1.0…2.0).
    pub p_exp: f32,
    /// Transport capacity slope exponent q (0.5…1.5).
    pub q_exp: f32,
    /// Bulk sediment density ρ (kg/m³).
    pub rho_sed: f32,
    /// Minimum slope used to avoid zero/unstable gradients.
    pub min_slope: f32,
    /// Diffusion substeps used to satisfy explicit stability (1..=10).
    pub subcycles: u32,
    /// If true, apply a flexural response immediately after deposition.
    pub couple_flexure: bool,
}

impl Default for SurfaceParams {
    fn default() -> Self {
        Self {
            // Soften stream-power coefficient to reduce aggressive incision
            k_stream: 1.0e-6,
            m_exp: 0.5,
            n_exp: 1.0,
            // Lower hillslope diffusivity for gentler smoothing
            k_diff: 0.05,
            k_tr: 0.1,
            p_exp: 1.3,
            q_exp: 1.0,
            rho_sed: 1800.0,
            // Slightly higher minimum slope to avoid near-zero unstable gradients
            min_slope: 2.0e-4,
            subcycles: 4,
            couple_flexure: false,
        }
    }
}

/// Diagnostic summary for one surface pass.
#[derive(Default, Clone, Copy, Debug)]
pub struct SurfaceStats {
    /// Total eroded volume (m³).
    pub eroded_m3: f64,
    /// Total deposited volume (m³).
    pub deposited_m3: f64,
    /// Volume residual (m³) = eroded − deposited.
    pub residual_m3: f64,
    /// Maximum erosion thickness at any cell (m).
    pub max_erosion_m: f32,
    /// Maximum deposition thickness at any cell (m).
    pub max_deposition_m: f32,
}

/// Apply surface processes in-place.
pub fn apply_surface_processes(
    grid: &Grid,
    _c: &[f32],
    depth_m: &mut [f32],
    sediment_m: &mut [f32],
    area_m2: &[f32],
    p: &SurfaceParams,
    dt_myr: f64,
) -> SurfaceStats {
    let n = grid.cells;
    if depth_m.len() < n || sediment_m.len() < n || area_m2.len() < n {
        return SurfaceStats::default();
    }
    if dt_myr <= 0.0 {
        return SurfaceStats::default();
    }

    let dt_yr = (dt_myr * 1.0e6) as f32;
    let dt_myr_f32 = dt_myr as f32;
    // Soft caps to prevent overly aggressive per-step changes
    const EROSION_CAP_M_PER_MYR: f32 = 0.5; // max erosion per Myr (m)
    const LAND_SUBMERGE_CAP_M_PER_MYR: f32 = 0.5; // max allowed land submergence per Myr (m)
    let erosion_cap_m = (EROSION_CAP_M_PER_MYR.max(0.0)) * dt_myr_f32;
    let land_submerge_cap_m = (LAND_SUBMERGE_CAP_M_PER_MYR.max(0.0)) * dt_myr_f32;

    // Elevation field (m)
    let depth_prev = depth_m.to_vec();
    let mut elev_m: Vec<f32> = depth_m.iter().map(|&d| -d).collect();
    // Priority-flood-like depression fill (iterative, land-only). This raises local closed pits
    // to their lowest spill elevation using repeated local minima elimination.
    // Ocean cells (elev<0) act as fixed outlets and are not modified.
    // More passes to reduce isolated pits that cause spurious steep local slopes
    fill_depressions_iterative(grid, &mut elev_m, 12);

    // 1) Local slope magnitude using least-squares plane fit in local EN frame.
    let slope = compute_slope_mag(grid, &elev_m);

    // 2) D8 steepest-descent receiver per cell; -1 indicates pit/outlet.
    let receiver: Vec<i32> = compute_steepest_receiver(grid, &elev_m);

    // 3) Flow accumulation A (m^2): start with own area then accumulate downstream.
    let mut area_acc: Vec<f32> = area_m2.to_vec();
    {
        let mut order: Vec<usize> = (0..n).collect();
        // Sort by descending elevation so upstream contributes before downstream
        order.sort_unstable_by(|&a, &b| {
            elev_m[b].partial_cmp(&elev_m[a]).unwrap_or(std::cmp::Ordering::Equal)
        });
        for i in order {
            let j = receiver[i];
            if j >= 0 {
                let jj = j as usize;
                area_acc[jj] += area_acc[i];
            }
        }
    }

    // 4) Fluvial erosion (only on land cells for MVP: depth<=0 → elev>=0)
    let mut erosion_thickness_m = vec![0.0f32; n];
    let mexp = p.m_exp.max(0.0);
    let nexp = p.n_exp.max(0.0);
    for i in 0..n {
        if depth_m[i] <= 0.0 {
            let s = slope[i].max(p.min_slope);
            let a = area_acc[i].max(1.0);
            let rate_m_per_yr = p.k_stream * a.powf(mexp) * s.powf(nexp);
            let de_raw = (rate_m_per_yr * dt_yr).max(0.0);
            // Cap per-step erosion thickness to avoid aggressive drawdown
            let de = if erosion_cap_m > 0.0 { de_raw.min(erosion_cap_m) } else { de_raw };
            erosion_thickness_m[i] = de;
        }
    }

    // Track eroded mass routed downstream (kg)
    let mut sed_mass_in = vec![0.0f64; n];
    let mut sed_mass_out = vec![0.0f64; n];
    let rho = p.rho_sed.max(1.0);
    let mut stats = SurfaceStats::default();

    // Convert erosion thickness to mass and send downstream (single receiver routing)
    {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| {
            elev_m[b].partial_cmp(&elev_m[a]).unwrap_or(std::cmp::Ordering::Equal)
        });
        for i in order {
            let thk = erosion_thickness_m[i];
            if thk > 0.0 {
                // Apply erosion locally to depth (Δdepth = +thk)
                depth_m[i] += thk;
                // Preferentially erode from prior sediment thickness if any
                let take_from_sed = thk.min(sediment_m[i]);
                sediment_m[i] -= take_from_sed;
                // Mass produced = volume * ρ
                let vol_m3 = (thk as f64) * (area_m2[i] as f64);
                stats.eroded_m3 += vol_m3;
                let mass_kg = vol_m3 * (rho as f64);
                let j = receiver[i];
                if j >= 0 {
                    sed_mass_in[j as usize] += mass_kg;
                } else {
                    // No receiver: accumulate locally for deposition test
                    sed_mass_in[i] += mass_kg;
                }
                stats.max_erosion_m = stats.max_erosion_m.max(thk);
            }
        }
    }

    // 5) Sediment transport capacity and deposition; propagate remaining mass downstream
    {
        let mut order: Vec<usize> = (0..n).collect();
        // Process from high to low so upstream supply is available first
        order.sort_unstable_by(|&a, &b| {
            elev_m[b].partial_cmp(&elev_m[a]).unwrap_or(std::cmp::Ordering::Equal)
        });
        let pexp = p.p_exp.max(0.0);
        let qexp = p.q_exp.max(0.0);
        for i in order {
            let s = slope[i].max(p.min_slope);
            let a = area_acc[i].max(1.0);
            // Capacity as mass per step (kg): (K_tr * A^p * S^q) * dt * ρ
            let cap_mass =
                ((p.k_tr * a.powf(pexp) * s.powf(qexp)) * dt_yr).max(0.0) as f64 * (rho as f64);
            let incoming = sed_mass_in[i];
            let carry = incoming.min(cap_mass);
            let deposit_mass = (incoming - carry).max(0.0);
            if deposit_mass > 0.0 {
                let dep_thk = (deposit_mass / (rho as f64)) / (area_m2[i] as f64);
                let dep_thk_f = dep_thk as f32;
                // Deposition raises surface → decreases depth
                depth_m[i] -= dep_thk_f;
                sediment_m[i] += dep_thk_f;
                stats.deposited_m3 += deposit_mass / (rho as f64);
                stats.max_deposition_m = stats.max_deposition_m.max(dep_thk_f);
            }
            // Route remaining mass downstream
            let j = receiver[i];
            if j >= 0 {
                sed_mass_in[j as usize] += carry;
            } else {
                // No receiver: deposit all remaining locally
                if carry > 0.0 {
                    let dep_thk = (carry / (rho as f64)) / (area_m2[i] as f64);
                    let dep_thk_f = dep_thk as f32;
                    depth_m[i] -= dep_thk_f;
                    sediment_m[i] += dep_thk_f;
                    stats.deposited_m3 += carry / (rho as f64);
                    stats.max_deposition_m = stats.max_deposition_m.max(dep_thk_f);
                }
            }
            sed_mass_out[i] = carry;
        }
    }

    // 6) Hillslope diffusion (resolution-aware): explicit update with distance-weighted Laplacian.
    // Stability (CFL): κ * dt * sum_j (1/dx_ij^2) <= 0.24 (safety vs 0.25).
    let mut smax = 0.0f64;
    {
        for i in 0..n {
            let lat_i = grid.latlon[i][0] as f64;
            let lon_i = grid.latlon[i][1] as f64;
            let cos_lat = lat_i.cos();
            let mut s = 0.0f64;
            for &nj in &grid.n1[i] {
                let j = nj as usize;
                let lat_j = grid.latlon[j][0] as f64;
                let lon_j = grid.latlon[j][1] as f64;
                let dx_e = (lon_j - lon_i) * cos_lat * 6_371_000.0f64;
                let dx_n = (lat_j - lat_i) * 6_371_000.0f64;
                let dx2 = (dx_e * dx_e + dx_n * dx_n).max(1.0);
                s += 1.0 / dx2;
            }
            if s > smax {
                smax = s;
            }
        }
    }
    let kappa = p.k_diff as f64;
    let mut subcycles = p.subcycles.max(1);
    if smax > 0.0 {
        let req = (kappa * (dt_yr as f64) * smax / 0.24).ceil() as u32;
        if req > subcycles {
            subcycles = req.min(128);
        }
    }
    let dt_sub = dt_yr / (subcycles as f32);
    for _ in 0..subcycles {
        let before: Vec<f32> = depth_m.to_vec();
        for i in 0..n {
            let lat_i = grid.latlon[i][0] as f64;
            let lon_i = grid.latlon[i][1] as f64;
            let cos_lat = lat_i.cos();
            let mut acc = 0.0f64;
            for &nj in &grid.n1[i] {
                let j = nj as usize;
                let lat_j = grid.latlon[j][0] as f64;
                let lon_j = grid.latlon[j][1] as f64;
                let dx_e = (lon_j - lon_i) * cos_lat * 6_371_000.0f64;
                let dx_n = (lat_j - lat_i) * 6_371_000.0f64;
                let dx2 = (dx_e * dx_e + dx_n * dx_n).max(1.0);
                let inv_dx2 = 1.0 / dx2;
                acc += (before[j] as f64 - before[i] as f64) * inv_dx2;
            }
            let d = (p.k_diff as f64) * (dt_sub as f64) * acc;
            depth_m[i] = (before[i] as f64 + d) as f32;
        }
    }

    // Shoreline stability: if a cell was land at start, allow at most a shallow submergence per step
    if land_submerge_cap_m >= 0.0 {
        for i in 0..n {
            if depth_prev[i] <= 0.0 && depth_m[i] > 0.0 {
                depth_m[i] = depth_m[i].min(land_submerge_cap_m);
            }
        }
    }
    stats.residual_m3 = stats.eroded_m3 - stats.deposited_m3;
    stats
}

fn compute_slope_mag(grid: &Grid, elev_m: &[f32]) -> Vec<f32> {
    // Fit plane z(x,y) via least-squares over 1-ring in local EN coordinates; slope = sqrt(gx^2+gy^2)
    const R: f64 = 6_371_000.0;
    let n = grid.cells;
    let mut slope = vec![0.0f32; n];
    for i in 0..n {
        let lat_i = grid.latlon[i][0] as f64;
        let lon_i = grid.latlon[i][1] as f64;
        let zi = elev_m[i] as f64;
        let cos_lat = lat_i.cos();
        let mut s_xx = 0.0;
        let mut s_xy = 0.0;
        let mut s_yy = 0.0;
        let mut s_xz = 0.0;
        let mut s_yz = 0.0;
        let mut used = 0usize;
        for &nj in &grid.n1[i] {
            let j = nj as usize;
            let lat_j = grid.latlon[j][0] as f64;
            let lon_j = grid.latlon[j][1] as f64;
            let dx = (lon_j - lon_i) * cos_lat * R; // east (m)
            let dy = (lat_j - lat_i) * R; // north (m)
            let dz = (elev_m[j] as f64) - zi;
            s_xx += dx * dx;
            s_xy += dx * dy;
            s_yy += dy * dy;
            s_xz += dx * dz;
            s_yz += dy * dz;
            used += 1;
        }
        if used >= 2 {
            let det = s_xx * s_yy - s_xy * s_xy;
            if det.abs() > 1e-12 {
                let inv_xx = s_yy / det;
                let inv_xy = -s_xy / det;
                let inv_yy = s_xx / det;
                let gx = inv_xx * s_xz + inv_xy * s_yz;
                let gy = inv_xy * s_xz + inv_yy * s_yz;
                slope[i] = ((gx * gx + gy * gy).sqrt()) as f32;
            } else {
                slope[i] = 0.0;
            }
        } else {
            slope[i] = 0.0;
        }
    }
    slope
}

fn compute_steepest_receiver(grid: &Grid, elev_m: &[f32]) -> Vec<i32> {
    let n = grid.cells;
    let mut recv = vec![-1i32; n];
    for i in 0..n {
        let zi = elev_m[i];
        let mut best_j: i32 = -1;
        let mut best_dz: f32 = 0.0;
        for &nj in &grid.n1[i] {
            let j = nj as usize;
            let dz = zi - elev_m[j]; // positive if neighbor is lower
            if dz > best_dz {
                best_dz = dz;
                best_j = j as i32;
            }
        }
        recv[i] = best_j;
    }
    recv
}

/// Iteratively fill depressions on land (elev>=0) by raising local minima to the
/// minimum neighbor elevation. Limited passes for stability; ocean (elev<0) remains unchanged.
fn fill_depressions_iterative(grid: &Grid, elev_m: &mut [f32], max_passes: usize) {
    let n = grid.cells;
    for _ in 0..max_passes {
        let mut changed = false;
        for i in 0..n {
            let zi = elev_m[i];
            if zi < 0.0 {
                continue; // ocean seed
            }
            let mut min_nb = f32::INFINITY;
            for &nj in &grid.n1[i] {
                let j = nj as usize;
                let zj = elev_m[j];
                if zj < min_nb {
                    min_nb = zj;
                }
            }
            if min_nb.is_finite() && zi < min_nb {
                elev_m[i] = min_nb;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
}
