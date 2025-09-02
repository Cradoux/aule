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

/// Solve eta using exact hypsometry prefix sums (monotone, bracketless).
///
/// Safety rails:
/// - Sanitizes inputs: clamp depths to [-9000, 11000], discard non-finite area.
/// - Bounds eta to ±20 km; on failure returns `eta_prev`.
pub fn solve_eta_hypsometry(depth_m: &[f32], area_m2: &[f32], v_target: f64, eta_prev: f64) -> f64 {
    let n = depth_m.len().min(area_m2.len());
    if n == 0 {
        return eta_prev;
    }
    // Pair and sanitize
    let mut cells: Vec<(f32, f32)> = Vec::with_capacity(n);
    for i in 0..n {
        let mut d = depth_m[i];
        if !d.is_finite() {
            d = 0.0;
        }
        d = d.clamp(-9000.0, 11000.0);
        let mut a = area_m2[i];
        if !a.is_finite() || a <= 0.0 {
            a = 0.0;
        }
        cells.push((d, a));
    }
    // Sort by depth ascending
    cells.sort_by(|(d1, _), (d2, _)| d1.partial_cmp(d2).unwrap_or(std::cmp::Ordering::Equal));
    // Prefix sums
    let mut s = 0.0f64; // area
    let mut p = 0.0f64; // area*depth
    let mut s_prefix: Vec<f64> = Vec::with_capacity(n);
    let mut p_prefix: Vec<f64> = Vec::with_capacity(n);
    for (d, a) in &cells {
        s += *a as f64;
        p += (*a as f64) * (*d as f64);
        s_prefix.push(s);
        p_prefix.push(p);
    }
    if s <= 0.0 {
        return eta_prev;
    }
    // Binary search index such that V(depth_k) <= V_target
    let mut lo = 0usize;
    let mut hi = n - 1;
    while lo < hi {
        let mid = (lo + hi) / 2;
        let d_mid = cells[mid].0 as f64;
        let v_mid = d_mid * s_prefix[mid] - p_prefix[mid];
        if v_mid <= v_target {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let k = lo.saturating_sub(1);
    let sk = s_prefix[k].max(1.0);
    let pk = p_prefix[k];
    let eta = (v_target + pk) / sk;
    let eta_bounded = eta.clamp(-20_000.0, 20_000.0);
    if eta_bounded.is_finite() {
        eta_bounded
    } else {
        eta_prev
    }
}

/// Compute eta for a target ocean area fraction using hypsometry (exact, area-weighted).
/// Returns (eta, ocean_area_m2, ocean_volume_m3).
pub fn solve_eta_for_ocean_area_fraction_hypso(
    depth_m: &[f32],
    area_m2: &[f32],
    target_ocean_frac: f64,
) -> (f64, f64, f64) {
    let n = depth_m.len().min(area_m2.len());
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    // Pair and sanitize
    let mut cells: Vec<(f32, f32)> = Vec::with_capacity(n);
    for i in 0..n {
        let mut d = depth_m[i];
        if !d.is_finite() {
            d = 0.0;
        }
        d = d.clamp(-9000.0, 11000.0);
        let mut a = area_m2[i];
        if !a.is_finite() || a <= 0.0 {
            a = 0.0;
        }
        cells.push((d, a));
    }
    // Sort by depth ascending
    cells.sort_by(|(d1, _), (d2, _)| d1.partial_cmp(d2).unwrap_or(std::cmp::Ordering::Equal));
    // Prefix sums
    let mut s = 0.0f64; // area
    let mut p = 0.0f64; // area*depth
    let mut s_prefix: Vec<f64> = Vec::with_capacity(n);
    let mut p_prefix: Vec<f64> = Vec::with_capacity(n);
    for (d, a) in &cells {
        s += *a as f64;
        p += (*a as f64) * (*d as f64);
        s_prefix.push(s);
        p_prefix.push(p);
    }
    if s <= 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let total_area = s;
    let target_ocean_area = (target_ocean_frac.clamp(0.0, 1.0)) * total_area;
    // Ocean area = total_area - S[k] for any eta in [depth_k, depth_{k+1}).
    // Find k so that total_area - S[k] ≈ target_ocean_area → S[k] ≈ total_area - target.
    let target_s = (total_area - target_ocean_area).clamp(0.0, total_area);
    let mut k = 0usize;
    while k + 1 < n && s_prefix[k] < target_s {
        k += 1;
    }
    // eta at this threshold; clamp within neighbor depths
    let eta = cells[k].0 as f64;
    let area_suffix = total_area - s_prefix[k];
    let p_suffix = p_prefix[n - 1] - p_prefix[k];
    let vol = (p_suffix - eta * area_suffix).max(0.0);
    (eta, area_suffix, vol)
}

/// Compute the reference ocean area and volume using elevation definition.
///
/// Ocean cells are those with (depth + eta_ref) > 0 (positive down). Returns a [`SeaLevelRef`]
/// capturing the ocean area (m^2) and volume (m^3) consistent with rendering and the solver.
pub fn compute_ref(depth_m: &[f32], area_m2: &[f32], eta_ref_m: f32) -> SeaLevelRef {
    assert_eq!(depth_m.len(), area_m2.len());
    let mut vol = 0.0f64;
    let mut area = 0.0f64;
    let eta = eta_ref_m as f64;
    for i in 0..depth_m.len() {
        // Ocean if (depth - eta) > 0 ⇔ z = η − depth < 0
        let d = depth_m[i] as f64 - eta;
        let a = area_m2[i] as f64;
        if d > 0.0 && a > 0.0 {
            vol += d * a;
            area += a;
        }
    }
    SeaLevelRef { volume_m3: vol, ocean_area_m2: area }
}

/// Set the world's sea-level reference to the current ocean state and return it.
///
/// This is used to re-baseline the global sea-level target after persistent changes
/// to topography (e.g., adding continents). If there are no ocean cells, sets both
/// fields to 0 and logs a note.
pub fn rebaseline(world: &mut World, area_m2: &[f32]) -> SeaLevelRef {
    let r = compute_ref(&world.depth_m, area_m2, world.sea.eta_m);
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
    // Bracket such that g(lo) <= 0 <= g(hi), where g(off) = area_frac(off) - target
    let mut lo = -12_000.0f64;
    let mut hi = 12_000.0f64;
    let mut a_lo = area_frac(lo) - target;
    let mut a_hi = area_frac(hi) - target;
    // Expand adaptively if needed (monotonic in off)
    let mut step = 12_000.0f64;
    let mut iters = 0;
    while !(a_lo <= 0.0 && a_hi >= 0.0) && iters < 32 {
        if a_lo > 0.0 {
            lo -= step;
            a_lo = area_frac(lo) - target;
        }
        if a_hi < 0.0 {
            hi += step;
            a_hi = area_frac(hi) - target;
        }
        step *= 2.0;
        iters += 1;
    }
    // If still not bracketed (degenerate cases), return the closer endpoint
    if !(a_lo <= 0.0 && a_hi >= 0.0) {
        return if a_lo.abs() < a_hi.abs() { lo } else { hi };
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

/// Solve η directly (meters) so that ocean_volume(depth+η) ~= target_vol_m3 via robust bracketing.
/// Starts from an initial guess `eta0_m` and expands the bracket until it encloses the root.
pub fn solve_eta_for_volume(
    depth_m: &[f32],
    area_m2: &[f32],
    target_vol_m3: f64,
    eta0_m: f64,
    tol_volume_m3: f64,
    max_iter: u32,
) -> f64 {
    #[inline]
    fn vol_eta(depth_m: &[f32], area_m2: &[f32], eta: f64) -> f64 {
        let mut v = 0.0f64;
        for (d, a) in depth_m.iter().zip(area_m2.iter()) {
            // Ocean thickness with η-convention (z = η − depth) is max(depth − η, 0)
            let dd = (*d as f64) - eta;
            if dd > 0.0 {
                v += dd * (*a as f64);
            }
        }
        v
    }

    // (hypsometry variant defined at top-level)

    // Center bracket around eta0
    let mut half_span = 6_000.0_f64;
    let mut lo = eta0_m - half_span;
    let mut hi = eta0_m + half_span;
    let mut f_lo = vol_eta(depth_m, area_m2, lo) - target_vol_m3;
    let mut f_hi = vol_eta(depth_m, area_m2, hi) - target_vol_m3;
    // Expand until we bracket or hit a generous limit
    let mut expands = 0;
    while !(f_lo <= 0.0 && f_hi >= 0.0) && expands < 16 {
        half_span *= 2.0;
        lo = eta0_m - half_span;
        hi = eta0_m + half_span;
        f_lo = vol_eta(depth_m, area_m2, lo) - target_vol_m3;
        f_hi = vol_eta(depth_m, area_m2, hi) - target_vol_m3;
        expands += 1;
    }
    // If still not bracketed, return the closer endpoint to avoid extreme jumps
    if !(f_lo <= 0.0 && f_hi >= 0.0) {
        return if f_lo.abs() < f_hi.abs() { lo } else { hi };
    }
    // Bisection
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..max_iter {
        mid = 0.5 * (lo + hi);
        let f_mid = vol_eta(depth_m, area_m2, mid) - target_vol_m3;
        if f_mid.abs() <= tol_volume_m3 || (hi - lo).abs() < 1e-6 {
            return mid;
        }
        if (f_mid > 0.0) == (f_lo > 0.0) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
    }
    mid
}

/// Solve uniform offset (meters) so that land fraction equals `target_land_frac`.
///
/// Land is defined where elevation > 0, with elevation = -(depth + off).
/// Positive `off` makes water deeper (more ocean); negative `off` raises land.
pub fn solve_offset_for_land_fraction(
    depth_m: &[f32],
    area_m2: &[f32],
    target_land_frac: f32,
    max_iter: u32,
) -> f64 {
    assert_eq!(depth_m.len(), area_m2.len());
    let total_area: f64 = area_m2.iter().map(|&a| a as f64).sum();
    if total_area <= 0.0 {
        return 0.0;
    }
    let target = (target_land_frac.clamp(0.0, 1.0)) as f64;
    let land_frac = |off: f64| -> f64 {
        let mut land = 0.0f64;
        for (&d, &a) in depth_m.iter().zip(area_m2.iter()) {
            if (-(d as f64 + off)) > 0.0 {
                land += a as f64;
            }
        }
        land / total_area
    };
    // Bracket off so that f(lo) >= target >= f(hi) (monotone decreasing in off)
    let mut lo = -8000.0f64;
    let mut hi = 8000.0f64;
    let mut f_lo = land_frac(lo);
    let mut f_hi = land_frac(hi);
    let mut step = 8000.0f64;
    let mut it = 0;
    while !(f_lo >= target && f_hi <= target) && it < 32 {
        if f_lo < target {
            lo -= step;
            f_lo = land_frac(lo);
        }
        if f_hi > target {
            hi += step;
            f_hi = land_frac(hi);
        }
        step *= 2.0;
        it += 1;
    }
    if !(f_lo >= target && f_hi <= target) {
        return if (f_lo - target).abs() < (f_hi - target).abs() { lo } else { hi };
    }
    // Bisection
    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        let f_mid = land_frac(mid);
        if (f_mid - target).abs() <= 1e-4 {
            // fraction tolerance
            return mid;
        }
        if f_mid > target {
            // too much land → move toward deeper water
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
