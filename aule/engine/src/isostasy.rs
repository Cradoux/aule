//! Global sea-level constraint utilities (constant ocean volume via uniform offset).

use crate::world::{SeaLevelRef, World};

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

/// Compute the reference ocean area and volume from the current depths.
///
/// Ocean cells are those with depth > 0 (positive down). Returns a [`SeaLevelRef`]
/// capturing the ocean area (m^2) and volume (m^3).
pub fn compute_ref(depth_m: &[f32], area_m2: &[f32]) -> SeaLevelRef {
    let (v, a) = ocean_volume_from_depth(depth_m, area_m2);
    SeaLevelRef { volume_m3: v, ocean_area_m2: a }
}

/// Set the world's sea-level reference to the current ocean state and return it.
///
/// This is used to re-baseline the global sea-level target after persistent changes
/// to topography (e.g., adding continents). If there are no ocean cells, sets both
/// fields to 0 and logs a note.
pub fn rebaseline(world: &mut World, area_m2: &[f32]) -> SeaLevelRef {
    let r = compute_ref(&world.depth_m, area_m2);
    world.sea_level_ref = Some(r);
    let area_sum: f64 = area_m2.iter().map(|&a| a as f64).sum();
    let frac = if area_sum > 0.0 { r.ocean_area_m2 / area_sum } else { 0.0 };
    if r.ocean_area_m2 == 0.0 {
        println!("[isostasy] rebaseline: no ocean cells; L disabled (offset=0)");
    } else {
        println!(
            "[isostasy] rebaseline: ocean_frac_ref={:.3} volume_ref={:.3e} m^3",
            frac, r.volume_m3
        );
    }
    r
}

/// In-place add a uniform sea-level offset (meters). Positive increases depth (down).
pub fn apply_sea_level_offset(depth_m: &mut [f32], offset_m: f64) {
    let off = offset_m as f32;
    for d in depth_m {
        *d += off;
    }
}

/// Bisection on uniform offset (meters) so ocean-area fraction == target_ocean_frac.
///
/// Convention: ocean if (depth_m[i] + off) > 0.0. Returns the offset in meters.
pub fn solve_offset_for_ocean_area_fraction(
    depth_m: &[f32],
    area_m2: &[f32],
    target_ocean_frac: f32,
    tol_frac: f64,
    max_iter: u32,
) -> f64 {
    assert_eq!(depth_m.len(), area_m2.len());
    let total_area: f64 = area_m2.iter().map(|&a| a as f64).sum();
    if total_area <= 0.0 {
        return 0.0;
    }
    let target = target_ocean_frac.clamp(0.0, 1.0) as f64;
    let area_frac = |off: f64| -> f64 {
        let mut ocean_area = 0.0f64;
        for (d, a) in depth_m.iter().zip(area_m2.iter()) {
            if (*d as f64 + off) > 0.0 {
                ocean_area += *a as f64;
            }
        }
        ocean_area / total_area
    };
    // Fixed bracket sufficient for global range
    let mut lo = -12_000.0f64;
    let mut hi = 12_000.0f64;
    let mut a_lo = area_frac(lo) - target;
    let mut a_hi = area_frac(hi) - target;
    // If bracket doesn't straddle, expand a few times conservatively
    let mut expand = 0;
    while a_lo > 0.0 && expand < 4 {
        lo -= 12_000.0;
        a_lo = area_frac(lo) - target;
        expand += 1;
    }
    expand = 0;
    while a_hi < 0.0 && expand < 4 {
        hi += 12_000.0;
        a_hi = area_frac(hi) - target;
        expand += 1;
    }
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..max_iter {
        mid = 0.5 * (lo + hi);
        let f_mid = area_frac(mid) - target;
        if f_mid.abs() <= tol_frac || (hi - lo).abs() < 1e-6 {
            return mid;
        }
        if (f_mid > 0.0) == (a_lo > 0.0) {
            lo = mid;
            a_lo = f_mid;
        } else {
            hi = mid;
        }
    }
    mid
}

/// Outer bisection on uplift amplitude (meters) over a continent template (0..1),
/// inner solve on sea-level offset to hit target land fraction.
/// Returns (amp_m, off_m).
#[allow(clippy::too_many_arguments)]
pub fn solve_amplitude_for_land_fraction(
    tpl: &[f32],
    base_depth_m: &[f32],
    area_m2: &[f32],
    target_land_frac: f32,
    amp_lo_m: f64,
    amp_hi_m: f64,
    tol_land: f64,
    tol_ocean: f64,
    max_iter: u32,
) -> (f64, f64) {
    assert_eq!(tpl.len(), base_depth_m.len());
    assert_eq!(tpl.len(), area_m2.len());
    let total_area: f64 = area_m2.iter().map(|&a| a as f64).sum();
    if total_area <= 0.0 {
        return (0.0, 0.0);
    }
    let target_land = target_land_frac.clamp(0.0, 1.0) as f64;
    let target_ocean = 1.0 - target_land;
    // Degenerate template → only sea-level offset
    if tpl.iter().all(|&t| t == 0.0) {
        let off = solve_offset_for_ocean_area_fraction(
            base_depth_m,
            area_m2,
            target_ocean as f32,
            tol_ocean,
            64,
        );
        return (0.0, off);
    }
    let mut lo = amp_lo_m.max(0.0);
    let mut hi = amp_hi_m.max(lo);
    let mut best = (lo, 0.0, 1.0f64); // (amp, off, land_frac)
    let mut tmp: Vec<f32> = vec![0.0; base_depth_m.len()];
    for _ in 0..max_iter {
        let amp = 0.5 * (lo + hi);
        // tmp = base - amp * tpl
        for ((out, &d), &t) in tmp.iter_mut().zip(base_depth_m.iter()).zip(tpl.iter()) {
            *out = d - (amp as f32) * t;
        }
        let off =
            solve_offset_for_ocean_area_fraction(&tmp, area_m2, target_ocean as f32, tol_ocean, 64);
        // Compute achieved land fraction with applied offset
        let mut ocean_area = 0.0f64;
        for (&d, &a) in tmp.iter().zip(area_m2.iter()) {
            if (d as f64 + off) > 0.0 {
                ocean_area += a as f64;
            }
        }
        let land = 1.0 - (ocean_area / total_area);
        best = (amp, off, land);
        let err = land - target_land;
        if err.abs() <= tol_land {
            return (amp, off);
        }
        if err > 0.0 {
            hi = amp;
        } else {
            lo = amp;
        }
    }
    (best.0, best.1)
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
