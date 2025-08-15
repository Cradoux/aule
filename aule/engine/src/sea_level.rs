//! Long-term eustatic sea-level controller (slow global offset on top of isostasy).

/// Eustasy policy that determines the additional global offset η(t).
#[derive(Clone, Copy, Debug)]
pub enum EustasyPolicy {
    /// Constant extra offset in meters (default 0)
    Constant {
        /// Constant η in meters applied each step
        eta_m: f64,
    },
    /// Track a target ocean fraction φ* (0..1)
    TargetOceanFraction {
        /// Target ocean area fraction φ* (0..1)
        phi_target: f64,
    },
    /// Linear trend in η over time
    LinearTrend {
        /// dη/dt in millimeters per year (mm/yr)
        d_eta_mm_per_yr: f64,
    },
    /// Cyclic variation of η
    Cyclic {
        /// Amplitude A of η (meters)
        amp_m: f64,
        /// Period P (Myr)
        period_myr: f64,
        /// Phase in degrees
        phase_deg: f64,
    },
}

/// Compute the ocean area fraction for a given uniform offset (depth > 0 is ocean).
fn ocean_fraction_with_offset(depth_m: &[f32], area_m2: &[f32], offset_m: f64) -> f64 {
    let mut ocean_area: f64 = 0.0;
    let mut total_area: f64 = 0.0;
    for (d, &a) in depth_m.iter().zip(area_m2.iter()) {
        let dd = (*d as f64) + offset_m;
        if dd > 0.0 {
            ocean_area += a as f64;
        }
        total_area += a as f64;
    }
    if total_area > 0.0 {
        ocean_area / total_area
    } else {
        0.0
    }
}

/// Solve for the uniform offset that achieves a target ocean fraction φ*.
/// Uses bisection on a wide bracket; tolerances are in meters of offset.
pub fn solve_offset_for_target_fraction(
    depth_m: &[f32],
    area_m2: &[f32],
    phi_target: f64,
    tol_m: f64,
    max_iter: usize,
) -> f64 {
    let mut lo = -10_000.0f64; // meters
    let mut hi = 10_000.0f64;
    let mut flo = ocean_fraction_with_offset(depth_m, area_m2, lo) - phi_target;
    let fhi0 = ocean_fraction_with_offset(depth_m, area_m2, hi) - phi_target;
    // If not bracketing due to extreme targets, just return the closer bound
    if flo.signum() == fhi0.signum() {
        return if flo.abs() <= fhi0.abs() { lo } else { hi };
    }
    let mut mid = 0.0f64;
    for _ in 0..max_iter {
        mid = 0.5 * (lo + hi);
        let fmid = ocean_fraction_with_offset(depth_m, area_m2, mid) - phi_target;
        if (hi - lo).abs() <= tol_m || fmid.abs() <= 1e-12 {
            break;
        }
        if fmid.signum() == flo.signum() {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
        }
    }
    mid
}

/// Update and return the eustatic offset η(t) according to the policy.
/// `offset_isostasy_m` is provided for context (already computed this step),
/// but policies do not depend on it in this MVP.
pub fn update_eustasy_eta(
    world_t_myr: f64,
    _dt_myr: f64,
    policy: &EustasyPolicy,
    depth_m: &[f32],
    area_m2: &[f32],
    _offset_isostasy_m: f64,
) -> f64 {
    match *policy {
        EustasyPolicy::Constant { eta_m } => eta_m,
        EustasyPolicy::TargetOceanFraction { phi_target } => {
            solve_offset_for_target_fraction(depth_m, area_m2, phi_target, 1e-3, 64)
        }
        EustasyPolicy::LinearTrend { d_eta_mm_per_yr } => {
            let d_per_yr_m = d_eta_mm_per_yr * 1.0e-3;
            world_t_myr * d_per_yr_m * 1.0e6
        }
        EustasyPolicy::Cyclic { amp_m, period_myr, phase_deg } => {
            if period_myr <= 0.0 {
                return 0.0;
            }
            let theta =
                2.0 * std::f64::consts::PI * (world_t_myr / period_myr) + phase_deg.to_radians();
            amp_m * theta.sin()
        }
    }
}

/// Returns offset `off_m` so that the ocean area fraction equals `target_ocean_frac`.
/// `depth_m`: positive is water depth (ocean); negative is land elevation.
pub fn solve_offset_for_ocean_area_fraction(
    depth_m: &[f32],
    area_m2: &[f32],
    target_ocean_frac: f32,
) -> f64 {
    let total: f64 = area_m2.iter().map(|a| *a as f64).sum();
    if total <= 0.0 {
        return 0.0;
    }
    let mut lo = -12_000.0_f64;
    let mut hi = 12_000.0_f64;
    let area = |off: f64| -> f64 {
        let mut acc = 0.0;
        for (d, a) in depth_m.iter().zip(area_m2.iter()) {
            if (*d as f64 + off) > 0.0 {
                acc += *a as f64;
            }
        }
        acc / total
    };
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let f = area(mid) - (target_ocean_frac as f64);
        if f.abs() < 1e-6 {
            return mid;
        }
        if f > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}
