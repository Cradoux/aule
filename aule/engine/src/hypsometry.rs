//! Area-weighted hypsometry histogram utilities.
//!
//! Inputs are per-cell depths (meters, positive down) and per-cell areas (m^2).
//! Returns bin centers (m) and area per bin (m^2), using `bins` uniformly-spaced bins
//! over `[min_m, max_m]` inclusive of the max edge.

/// Compute an area-weighted histogram of depths.
///
/// - `depth_m`: per-cell depths in meters (positive down)
/// - `area_m2`: per-cell areas in square meters
/// - `bins`: number of bins (>= 1)
/// - `min_m`, `max_m`: histogram domain (meters). If equal or invalid, returns empty outputs.
///
/// Deterministic and allocation-free per element; no unsafe.
pub fn histogram_area_weighted(
    depth_m: &[f32],
    area_m2: &[f32],
    bins: usize,
    min_m: f32,
    max_m: f32,
) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(depth_m.len(), area_m2.len());
    if bins == 0 {
        return (Vec::new(), Vec::new());
    }
    if max_m.partial_cmp(&min_m) != Some(std::cmp::Ordering::Greater) {
        return (Vec::new(), Vec::new());
    }
    let nb = bins;
    let mut area_per_bin: Vec<f64> = vec![0.0; nb];
    let width = (max_m - min_m) as f64 / (nb as f64);
    if width == 0.0 {
        return (Vec::new(), Vec::new());
    }
    // Accumulate area per bin with clamped indexing
    for i in 0..depth_m.len() {
        let d = depth_m[i] as f64;
        let a = area_m2[i] as f64;
        let t = ((d - min_m as f64) / (max_m as f64 - min_m as f64)) * (nb as f64);
        if t.is_nan() || !t.is_finite() {
            continue;
        }
        let mut idx = t.floor() as isize;
        if idx < 0 {
            idx = 0;
        }
        if idx >= nb as isize {
            idx = (nb as isize) - 1;
        }
        area_per_bin[idx as usize] += a;
    }
    // Centers
    let mut centers: Vec<f64> = Vec::with_capacity(nb);
    for b in 0..nb {
        let c = (min_m as f64) + (b as f64 + 0.5) * width;
        centers.push(c);
    }
    (centers, area_per_bin)
}
