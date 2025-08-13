//! Global sea-level constraint utilities (constant ocean volume via uniform offset).

/// Compute ocean volume (m^3) and ocean area (m^2) from depths and cell areas.
/// depth_m: positive down; only depths >= 0 contribute volume/area.
pub fn ocean_volume_from_depth(depth_m: &[f32], area_m2: &[f32]) -> (f64, f64) {
    assert_eq!(depth_m.len(), area_m2.len());
    let mut vol = 0.0f64;
    let mut area = 0.0f64;
    for i in 0..depth_m.len() {
        let d = depth_m[i] as f64;
        let a = area_m2[i] as f64;
        if d > 0.0 && a > 0.0 {
            vol += d * a;
            area += a;
        }
    }
    (vol, area)
}

/// In-place add a uniform sea-level offset (meters). Positive increases depth (down).
pub fn apply_sea_level_offset(depth_m: &mut [f32], offset_m: f64) {
    let off = offset_m as f32;
    for d in depth_m {
        *d += off;
    }
}

/// Solve for offset (m) so that ocean_volume(depth+offset) ~= target_volume_m3 via bisection.
pub fn solve_offset_for_volume(
    depth_m: &[f32],
    area_m2: &[f32],
    target_vol_m3: f64,
    tol_volume_m3: f64,
    max_iter: u32,
) -> f64 {
    #[inline]
    fn vol_off(depth_m: &[f32], area_m2: &[f32], off: f64) -> f64 {
        let mut v = 0.0f64;
        for (d, a) in depth_m.iter().zip(area_m2.iter()) {
            let dd = (*d as f64) + off;
            if dd > 0.0 {
                v += dd * (*a as f64);
            }
        }
        v
    }

    // Initial bracket
    let mut lo = -10_000.0_f64;
    let mut hi = 10_000.0_f64;
    let mut f_lo = vol_off(depth_m, area_m2, lo) - target_vol_m3;
    let mut _f_hi = vol_off(depth_m, area_m2, hi) - target_vol_m3;

    // Enlarge bracket if needed
    let mut expand = 0;
    while f_lo > 0.0 && expand < 4 {
        lo -= 10_000.0;
        f_lo = vol_off(depth_m, area_m2, lo) - target_vol_m3;
        expand += 1;
    }
    expand = 0;
    while _f_hi < 0.0 && expand < 4 {
        hi += 10_000.0;
        _f_hi = vol_off(depth_m, area_m2, hi) - target_vol_m3;
        expand += 1;
    }

    // Bisection
    let mut off = 0.5 * (lo + hi);
    for _ in 0..max_iter {
        off = 0.5 * (lo + hi);
        let f_mid = vol_off(depth_m, area_m2, off) - target_vol_m3;
        if f_mid.abs() <= tol_volume_m3 || (hi - lo).abs() < 1e-6 {
            return off;
        }
        if (f_mid > 0.0) == (f_lo > 0.0) {
            lo = off;
            f_lo = f_mid;
        } else {
            hi = off;
        }
    }
    off
}
