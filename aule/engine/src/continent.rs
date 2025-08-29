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

/// Simple description of a great-circle belt used for inherited orogeny imprints.
/// The belt is centered on the great circle defined by the plane normal `n_hat`
/// (unit), with a half-width specified in kilometers.
#[derive(Clone, Copy, Debug)]
pub struct GreatCircleBelt {
    /// Unit plane normal describing the great circle (|n_hat| == 1)
    pub n_hat: [f64; 3],
    /// Half-width of the belt in kilometers
    pub half_width_km: f64,
    /// Peak uplift magnitude (negative makes land higher) in meters
    pub uplift_m: f32,
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

/// Build a supercontinent-like template by placing a ribbon of caps along an
/// equatorial great circle with several lobes. This is deliberately simple but
/// produces a contiguous landmass suitable as a starting point.
pub fn build_supercontinent_template(
    grid: &Grid,
    seed: u64,
    n_lobes: u32,
    lobe_radius_km: f64,
    falloff_km: f64,
) -> Vec<f32> {
    let mut template: Vec<f32> = vec![0.0; grid.cells];
    // Deterministic RNG
    let ns: u64 = 0x7375706572636e74; // "supercnt"
    let mut rng = StdRng::seed_from_u64(seed ^ ns);
    // Choose an equatorial great circle rotated by a random yaw about Z and random tilt
    let yaw = rng.gen_range(0.0..(2.0 * std::f64::consts::PI));
    let tilt = rng.gen_range(-0.6..0.6); // radians
    let rot_yaw = rotation_z(yaw);
    let rot_tilt = rotation_x(tilt);
    let r_plate_m = 6_371_000.0f64;
    let sigma_rad = (falloff_km * 1000.0) / r_plate_m;
    let r0_rad = (lobe_radius_km * 1000.0) / r_plate_m;
    let lobes = n_lobes.max(1) as usize;
    // Evenly spaced lobes across +/- 120 degrees along the ribbon
    let span = 2.0 * std::f64::consts::PI * (120.0 / 360.0);
    let start = -0.5 * span;
    let step = if lobes > 1 { span / ((lobes - 1) as f64) } else { 1.0 };
    let mut centers: Vec<[f64; 3]> = Vec::with_capacity(lobes);
    for i in 0..lobes {
        let lon = start + (i as f64) * step;
        let lat = 0.0f64;
        let mut v = lon_lat_to_xyz(lon, lat);
        v = mat3_mul(rot_tilt, v);
        v = mat3_mul(rot_yaw, v);
        centers.push(v);
    }
    for (i, r) in grid.pos_xyz.iter().enumerate() {
        let rhat = [r[0] as f64, r[1] as f64, r[2] as f64];
        let mut inv_prod = 1.0f64;
        for c in &centers {
            let dot = (rhat[0] * c[0] + rhat[1] * c[1] + rhat[2] * c[2]).clamp(-1.0, 1.0);
            let theta = dot.acos();
            let tk = if theta <= r0_rad { 1.0 } else { (-0.5 * (theta / sigma_rad).powi(2)).exp() };
            inv_prod *= 1.0 - tk;
        }
        template[i] = (1.0 - inv_prod) as f32;
    }
    template
}

/// Imprint inherited mountain belts as shallow negative-depth bands along great circles.
/// Returns a per-cell delta that can be added to `depth_m` (negative is uplift).
pub fn imprint_orogenic_belts(grid: &Grid, belts: &[GreatCircleBelt]) -> Vec<f32> {
    let mut delta = vec![0.0f32; grid.cells];
    let r_m = 6_371_000.0f64;
    for (i, r) in grid.pos_xyz.iter().enumerate() {
        let rhat = [r[0] as f64, r[1] as f64, r[2] as f64];
        for b in belts {
            let n = normalize3(b.n_hat);
            // Angular distance to great circle plane = |asin(n·r)|
            let ang = (dot3(n, rhat)).abs().asin();
            let hw_rad = (b.half_width_km * 1000.0) / r_m;
            if ang <= hw_rad {
                // Smooth taper to edges using cosine bell
                let t = (ang / hw_rad).clamp(0.0, 1.0) as f32;
                let w = 0.5 * (1.0 + (std::f32::consts::PI * (t - 1.0)).cos());
                delta[i] += -b.uplift_m * w;
            }
        }
    }
    delta
}

/// Parameters for supercontinent inherited belts.
#[derive(Clone, Copy, Debug)]
pub struct BeltParams {
    /// Half-width of the primary belt (km)
    pub half_width_km_primary: f64,
    /// Half-width of the secondary belt (km)
    pub half_width_km_secondary: f64,
    /// Peak uplift magnitude of primary belt (m; negative used when applied to depth)
    pub uplift_primary_m: f32,
    /// Peak uplift magnitude of secondary belt (m; negative used when applied to depth)
    pub uplift_secondary_m: f32,
    /// Diagonal angle in degrees between primary and secondary belts in the continent frame
    pub diag_angle_deg: f64,
}

impl Default for BeltParams {
    fn default() -> Self {
        Self {
            half_width_km_primary: 350.0,
            half_width_km_secondary: 220.0,
            uplift_primary_m: 500.0,
            uplift_secondary_m: 300.0,
            diag_angle_deg: 35.0,
        }
    }
}

/// Build great-circle belts aligned to the supercontinent ribbon orientation used in
/// `build_supercontinent_template`, derived deterministically from `seed`.
pub fn build_supercontinent_belts(seed: u64, p: BeltParams) -> Vec<GreatCircleBelt> {
    // Recreate the same orientation RNG draws as the template to align belts
    let ns: u64 = 0x7375706572636e74; // "supercnt"
    let mut rng = StdRng::seed_from_u64(seed ^ ns);
    let yaw = rng.gen_range(0.0..(2.0 * std::f64::consts::PI));
    let tilt = rng.gen_range(-0.6..0.6); // radians
    let rot_yaw = rotation_z(yaw);
    let rot_tilt = rotation_x(tilt);
    // Transform basis axes into continent frame
    let ex = [1.0, 0.0, 0.0];
    let ey = [0.0, 1.0, 0.0];
    let mut exp = mat3_mul(rot_tilt, ex);
    exp = mat3_mul(rot_yaw, exp);
    let mut eyp = mat3_mul(rot_tilt, ey);
    eyp = mat3_mul(rot_yaw, eyp);
    let n1 = normalize3(exp);
    // Secondary belt is a diagonal combination within the continent frame
    let phi = (p.diag_angle_deg.to_radians()).clamp(-std::f64::consts::PI, std::f64::consts::PI);
    let n2u = [
        phi.cos() * n1[0] + phi.sin() * eyp[0],
        phi.cos() * n1[1] + phi.sin() * eyp[1],
        phi.cos() * n1[2] + phi.sin() * eyp[2],
    ];
    let n2 = normalize3(n2u);
    vec![
        GreatCircleBelt { n_hat: n1, half_width_km: p.half_width_km_primary, uplift_m: p.uplift_primary_m },
        GreatCircleBelt { n_hat: n2, half_width_km: p.half_width_km_secondary, uplift_m: p.uplift_secondary_m },
    ]
}

// ----------------------- small 3D helpers (double precision) -----------------------
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let mut s = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if s <= 0.0 {
        s = 1.0;
    }
    [v[0] / s, v[1] / s, v[2] / s]
}

#[inline]
fn lon_lat_to_xyz(lon: f64, lat: f64) -> [f64; 3] {
    let cl = lat.cos();
    [cl * lon.cos(), lat.sin(), -cl * lon.sin()]
}

#[inline]
fn rotation_z(theta: f64) -> [[f64; 3]; 3] {
    let (c, s) = (theta.cos(), theta.sin());
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

#[inline]
fn rotation_x(theta: f64) -> [[f64; 3]; 3] {
    let (c, s) = (theta.cos(), theta.sin());
    [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
}

#[inline]
fn mat3_mul(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
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

/// Advect the continental fraction `C` and thickness `th_c_m` using a simple
/// semi-Lagrangian nearest-neighbour backtrace on the lat/lon grid.
///
/// This is a lightweight, deterministic CPU MVP intended to make motion visible.
/// It backtraces from each cell center by `v * dt` on a spherical lat/lon chart
/// and samples from the nearest of the 1-ring neighbours (including self).
pub fn advect_c_thc(
    grid: &Grid,
    vel_en_m_per_yr: &[[f32; 2]],
    dt_myr: f64,
    c: &mut [f32],
    th_c_m: &mut [f32],
) {
    let n = grid.cells;
    if c.len() != n || th_c_m.len() != n || vel_en_m_per_yr.len() != n {
        return;
    }
    // Precompute backtraced targets in lat/lon
    const R_EARTH_M: f64 = 6_371_000.0;
    let years = dt_myr.max(0.0) * 1.0e6;
    let mut c_new = vec![0.0f32; n];
    let mut thc_new = vec![0.0f32; n];
    for i in 0..n {
        let lat = grid.latlon[i][0] as f64;
        let lon = grid.latlon[i][1] as f64;
        let ve = vel_en_m_per_yr[i][0] as f64;
        let vn = vel_en_m_per_yr[i][1] as f64;
        let dx = ve * years; // meters east
        let dy = vn * years; // meters north
        let dlat = dy / R_EARTH_M;
        let dlon = if lat.abs() < std::f64::consts::FRAC_PI_2 {
            dx / (R_EARTH_M * lat.cos().max(1e-9))
        } else {
            0.0
        };
        // Backtrace: source position = current - displacement
        let src_lat = lat - dlat;
        let mut src_lon = lon - dlon;
        // Wrap lon to [-pi,pi]
        if src_lon > std::f64::consts::PI {
            src_lon -= 2.0 * std::f64::consts::PI;
        }
        if src_lon < -std::f64::consts::PI {
            src_lon += 2.0 * std::f64::consts::PI;
        }

        // Choose nearest among self + 1-ring neighbours in lat/lon
        let mut best_j = i;
        let mut best_d2 = (grid.latlon[i][0] as f64 - src_lat).powi(2)
            + (grid.latlon[i][1] as f64 - src_lon).powi(2);
        for &nj in &grid.n1[i] {
            let j = nj as usize;
            let d2 = (grid.latlon[j][0] as f64 - src_lat).powi(2)
                + (grid.latlon[j][1] as f64 - src_lon).powi(2);
            if d2 < best_d2 {
                best_d2 = d2;
                best_j = j;
            }
        }
        c_new[i] = c[best_j];
        thc_new[i] = th_c_m[best_j];
    }
    c.copy_from_slice(&c_new);
    th_c_m.copy_from_slice(&thc_new);
}

/// Apply uplift to `depth_m` using the fields `C` (0..1) and `th_c_m` (m).
/// Negative uplift makes water shallower (land higher). This simply applies:
/// depth[i] += -(C[i] * th_c_m[i]).
pub fn apply_uplift_from_c_thc(depth_m: &mut [f32], c: &[f32], th_c_m: &[f32]) {
    let n = depth_m.len().min(c.len()).min(th_c_m.len());
    for i in 0..n {
        depth_m[i] += -(c[i] * th_c_m[i]);
    }
}

/// Build a cratonic thickness field (meters) from a continent template with tapered core and mild spatial noise.
///
/// Deterministic given `seed`. Thickness = base_km + extra_km * t^p + noise_km * t * n_smooth,
/// where t ∈ [0,1] is the template value and n_smooth is a single-ring averaged hash noise in [-1,1].
pub fn build_craton_thickness(
    grid: &Grid,
    template: &[f32],
    base_km: f32,
    extra_km: f32,
    noise_km: f32,
    seed: u64,
    power: f32,
) -> Vec<f32> {
    let n = grid.cells.min(template.len());
    // Hash-based per-cell noise in [-1,1]
    fn hash32(mut x: u32) -> u32 {
        x ^= x >> 16;
        x = x.wrapping_mul(0x7feb352d);
        x ^= x >> 15;
        x = x.wrapping_mul(0x846ca68b);
        x ^= x >> 16;
        x
    }
    let seed_lo = (seed as u32) ^ 0xC0A71C0E;
    let mut noise_raw: Vec<f32> = vec![0.0; n];
    for (i, nr) in noise_raw.iter_mut().enumerate().take(n) {
        let h = hash32(i as u32 ^ seed_lo);
        let u = (h as f32) * (1.0 / 4294967295.0); // [0,1]
        *nr = 2.0 * u - 1.0; // [-1,1]
    }
    // Single-ring smoothing to reduce speckle while keeping determinism and low cost
    let mut noise_s: Vec<f32> = vec![0.0; n];
    for i in 0..n {
        let mut sum = noise_raw[i];
        let mut cnt = 1.0f32;
        for &nj in &grid.n1[i] {
            let j = nj as usize;
            if j < n {
                sum += noise_raw[j];
                cnt += 1.0;
            }
        }
        noise_s[i] = sum / cnt;
    }
    // Map template to thickness (km), then to meters
    let mut th_m: Vec<f32> = vec![0.0; n];
    for i in 0..n {
        let t = template[i].clamp(0.0, 1.0);
        let core = t.powf(power.max(1.0));
        let th_km = base_km + extra_km * core + noise_km * t * noise_s[i];
        th_m[i] = (th_km.max(0.0)) * 1000.0;
    }
    th_m
}
