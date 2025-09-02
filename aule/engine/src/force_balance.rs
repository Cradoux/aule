//! Simple force-balance update for plate Euler pole angular speeds.
//!
//! This MVP aggregates boundary kinematics per plate and nudges each plate's
//! angular speed `omega_rad_yr` toward a balance between driving (slab/ridge/transform)
//! and a linear damping term.

use crate::{boundaries::Boundaries, grid::Grid, plates::Plates};

/// Parameters for force-balance Euler pole updates.
#[derive(Clone, Copy, Debug)]
pub struct FbParams {
    /// Gain factor converting integrated boundary drive to Δω per step [rad/yr per unit].
    pub gain: f32,
    /// Linear damping coefficient on omega [1/Myr] (scaled by dt_myr).
    pub damp_per_myr: f32,
    /// Weight for convergent normal speed (slab pull proxy).
    pub k_conv: f32,
    /// Weight for divergent normal opening (ridge push proxy).
    pub k_div: f32,
    /// Weight for transform tangential shear (fault drag proxy; negative drive if desired).
    pub k_trans: f32,
    /// Maximum change of omega per step (rad/yr).
    pub max_domega: f32,
    /// Clamp absolute omega to this maximum (rad/yr).
    pub max_omega: f32,
}

impl Default for FbParams {
    fn default() -> Self {
        Self {
            gain: 1.0e-12,
            damp_per_myr: 0.2,
            k_conv: 1.0,
            k_div: 0.5,
            k_trans: 0.1,
            max_domega: 5.0e-9,
            max_omega: 2.0e-7,
        }
    }
}

/// Apply a single force-balance Euler speed update in-place on `plates.omega_rad_yr`.
///
/// - Accumulates a boundary drive per plate by summing edge-class-weighted velocities
///   times an edge-length proxy.
/// - Updates each plate's `omega_rad_yr` via Δω = gain * drive − damp * ω, with clamps.
pub fn apply_force_balance(
    grid: &Grid,
    bounds: &Boundaries,
    plate_id: &[u16],
    plates: &mut Plates,
    area_m2: &[f32],
    dt_myr: f64,
    p: FbParams,
) {
    let nplates = plates.pole_axis.len();
    if nplates == 0 {
        return;
    }
    let mut drive = vec![0.0f64; nplates];

    // Approximate edge length scale from cell areas.
    let mut cell_len = vec![0.0f64; grid.cells];
    for i in 0..grid.cells {
        cell_len[i] = (area_m2[i] as f64).sqrt();
    }

    for ek in &bounds.edge_kin {
        let u = ek.u as usize;
        let v = ek.v as usize;
        let pu = plate_id[u] as usize;
        let pv = plate_id[v] as usize;
        if pu == pv {
            continue;
        }
        let w = 0.5 * (cell_len[u] + cell_len[v]);
        match ek.class as u8 {
            2 => {
                // convergent: use closing speed magnitude
                let fnorm = (-ek.n_m_per_yr).max(0.0) as f64;
                let d = (p.k_conv as f64) * fnorm * w;
                drive[pu] += d;
                drive[pv] += d;
            }
            1 => {
                // divergent: opening
                let fnorm = (ek.n_m_per_yr).max(0.0) as f64;
                let d = (p.k_div as f64) * fnorm * w;
                drive[pu] += d;
                drive[pv] += d;
            }
            3 => {
                // transform: shear magnitude
                let ts = (ek.t_m_per_yr.abs()) as f64;
                let d = (p.k_trans as f64) * ts * w;
                drive[pu] += d;
                drive[pv] += d;
            }
            _ => {}
        }
    }

    let dt = dt_myr.max(0.0);
    let damp = (p.damp_per_myr as f64) * dt;
    let gain = (p.gain as f64) * dt;
    let max_do = (p.max_domega as f64).max(0.0);
    let max_o = (p.max_omega as f64).max(0.0);
    for (pid, _) in drive.iter().enumerate().take(nplates) {
        let omega = plates.omega_rad_yr[pid] as f64;
        // d(omega) = gain * drive - damp * omega
        let mut domega = gain * drive[pid] - damp * omega;
        if domega > max_do {
            domega = max_do;
        }
        if domega < -max_do {
            domega = -max_do;
        }
        let mut new_o = omega + domega;
        if new_o < 0.0 {
            new_o = 0.0;
        }
        if new_o > max_o {
            new_o = max_o;
        }
        plates.omega_rad_yr[pid] = new_o as f32;
    }
}
