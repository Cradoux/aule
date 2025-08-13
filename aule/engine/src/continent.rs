//! Synthetic continental mask (few large cratons) built procedurally on the unit sphere.
//! Deterministic given seed; uses smooth-union of Gaussian caps with a flat plateau.

use crate::grid::Grid;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Parameters controlling continent generation and application.
#[derive(Clone, Copy)]
pub struct ContinentParams {
    /// RNG seed (deterministic)
    pub seed: u64,
    /// Number of continents (caps)
    pub n_continents: u32,
    /// Plateau mean radius (km)
    pub mean_radius_km: f64,
    /// Gaussian falloff sigma (km)
    pub falloff_km: f64,
    /// Plateau uplift (negative is up/shallower)
    pub plateau_uplift_m: f32,
    /// Optional target land fraction (0..1); None → manual amplitude
    pub target_land_fraction: Option<f64>,
}

/// Outputs for continents
pub struct ContinentField {
    /// Unitless 0..1 template per cell; will be scaled by amplitude
    pub uplift_template_m: Vec<f32>,
    /// Applied uplift in meters (negative up) per cell
    pub uplift_applied_m: Vec<f32>,
    /// Land mask after application (depth ≤ 0)
    pub mask_land: Vec<bool>,
    /// Area-weighted fraction of land
    pub land_fraction: f64,
}

/// Build continent template (unitless 0..1) from smooth union of caps.
pub fn build_continents(grid: &Grid, p: ContinentParams) -> ContinentField {
    let mut template: Vec<f32> = vec![0.0; grid.cells];

    // Deterministic RNG (namespaced)
    let ns: u64 = 0x636f6e74696e65; // "contine"
    let mut rng = StdRng::seed_from_u64(p.seed ^ ns);

    // Generate cap centers uniformly on sphere
    let mut centers: Vec<[f64; 3]> = Vec::with_capacity(p.n_continents as usize);
    for _ in 0..p.n_continents {
        // Uniform on sphere via z in [-1,1], phi in [0,2pi)
        let z: f64 = rng.gen_range(-1.0..=1.0);
        let phi: f64 = rng.gen_range(0.0..(2.0 * std::f64::consts::PI));
        let rxy = (1.0 - z * z).sqrt();
        centers.push([rxy * phi.cos(), rxy * phi.sin(), z]);
    }

    // Parameters in radians
    const R_EARTH_M: f64 = 6_371_000.0;
    let sigma_rad = (p.falloff_km * 1000.0) / R_EARTH_M;
    let r0_rad = (p.mean_radius_km * 1000.0) / R_EARTH_M;

    for (i, r) in grid.pos_xyz.iter().enumerate() {
        let rhat = [r[0] as f64, r[1] as f64, r[2] as f64];
        let mut inv_prod = 1.0f64; // ∏(1 - t_k)
        for c in &centers {
            let dot = rhat[0] * c[0] + rhat[1] * c[1] + rhat[2] * c[2];
            let dot = dot.clamp(-1.0, 1.0);
            let theta = dot.acos();
            let tk = if theta <= r0_rad {
                1.0
            } else {
                let x = theta / sigma_rad;
                (-0.5 * x * x).exp()
            };
            inv_prod *= 1.0 - tk;
        }
        let t = 1.0 - inv_prod;
        template[i] = t as f32;
    }

    ContinentField {
        uplift_template_m: template,
        uplift_applied_m: vec![0.0; grid.cells],
        mask_land: vec![false; grid.cells],
        land_fraction: 0.0,
    }
}

/// Solve amplitude (meters, ≥0) to reach target land fraction via bisection.
pub fn solve_amplitude_for_target_land_fraction(
    depth_before_m: &[f32],
    uplift_template_m: &[f32],
    area_m2: &[f32],
    target_frac: f64,
    tol: f64,
    max_iter: u32,
) -> f32 {
    let total_area: f64 = area_m2.iter().map(|&a| a as f64).sum();
    if total_area <= 0.0 {
        return 0.0;
    }
    let mut lo = 0.0f64;
    let mut hi = 10_000.0f64; // 10 km uplift

    let f = |amp: f64| -> f64 {
        let mut land_area = 0.0f64;
        for i in 0..depth_before_m.len() {
            let d = depth_before_m[i] as f64 - amp * (uplift_template_m[i] as f64);
            if d <= 0.0 {
                land_area += area_m2[i] as f64;
            }
        }
        (land_area / total_area) - target_frac
    };

    // Ensure bracket contains root
    let mut f_lo = f(lo);
    let mut _f_hi = f(hi);
    let mut expand = 0;
    while f_lo > 0.0 && expand < 4 {
        hi *= 2.0;
        _f_hi = f(hi);
        expand += 1;
    }
    expand = 0;
    while _f_hi < 0.0 && expand < 4 {
        hi *= 2.0;
        _f_hi = f(hi);
        expand += 1;
    }

    let mut amp = 0.5 * (lo + hi);
    for _ in 0..max_iter {
        amp = 0.5 * (lo + hi);
        let fm = f(amp);
        if fm.abs() <= tol || (hi - lo).abs() < 1e-6 {
            return amp as f32;
        }
        if (fm > 0.0) == (f_lo > 0.0) {
            lo = amp;
            f_lo = fm;
        } else {
            hi = amp;
        }
    }
    amp as f32
}

/// Apply continents uplift: depth += -(amp * template). Returns (mask_land, land_fraction).
pub fn apply_continents(
    depth_m: &mut [f32],
    uplift_template_m: &[f32],
    amp_m: f32,
    area_m2: &[f32],
) -> (Vec<bool>, f64) {
    for i in 0..depth_m.len() {
        depth_m[i] += -(amp_m * uplift_template_m[i]);
    }
    let total_area: f64 = area_m2.iter().map(|&a| a as f64).sum();
    let mut land_area = 0.0f64;
    let mut mask = vec![false; depth_m.len()];
    for i in 0..depth_m.len() {
        if depth_m[i] <= 0.0 {
            land_area += area_m2[i] as f64;
            mask[i] = true;
        }
    }
    let frac = if total_area > 0.0 { land_area / total_area } else { 0.0 };
    (mask, frac)
}
