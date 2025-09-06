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
    let pc = crate::PhysConsts::default();
    let r_earth_m = pc.r_earth_m;
    let sigma_rad = (p.falloff_km * 1000.0) / r_earth_m;
    let r0_rad = (p.mean_radius_km * 1000.0) / r_earth_m;

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
        GreatCircleBelt {
            n_hat: n1,
            half_width_km: p.half_width_km_primary,
            uplift_m: p.uplift_primary_m,
        },
        GreatCircleBelt {
            n_hat: n2,
            half_width_km: p.half_width_km_secondary,
            uplift_m: p.uplift_secondary_m,
        },
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
    let pc = crate::PhysConsts::default();
    let r_earth_m = pc.r_earth_m;
    let years = dt_myr.max(0.0) * 1.0e6;
    let mut c_new = vec![0.0f32; n];
    let mut thc_new = vec![0.0f32; n];
    // Aggregate diagnostics instead of per-cell spam
    let mut warn_count: u32 = 0;
    let mut warn_min: f32 = f32::MAX;
    let mut warn_max: f32 = f32::MIN;
    let mut warn_sum: f64 = 0.0;
    for i in 0..n {
        let lat = grid.latlon[i][0] as f64;
        let lon = grid.latlon[i][1] as f64;
        let ve = vel_en_m_per_yr[i][0] as f64;
        let vn = vel_en_m_per_yr[i][1] as f64;
        let dx = ve * years; // meters east
        let dy = vn * years; // meters north
        let dlat = dy / r_earth_m;
        let dlon = if lat.abs() < std::f64::consts::FRAC_PI_2 {
            dx / (r_earth_m * lat.cos().max(1e-9))
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

        // Conservative interpolation over 1-ring + 2-ring: weights ~ 1/d2, normalized
        let mut wsum = 1.0f64; // include self with unit weight
        let mut c_acc = c[i] as f64;
        let mut thc_acc = th_c_m[i] as f64;
        let add_neighbor = |j: usize, wsum: &mut f64, c_acc: &mut f64, thc_acc: &mut f64| {
            let d2 = (grid.latlon[j][0] as f64 - src_lat).powi(2)
                + (grid.latlon[j][1] as f64 - src_lon).powi(2);
            let w = 1.0 / d2.max(1e-12);
            *wsum += w;
            *c_acc += w * (c[j] as f64);
            *thc_acc += w * (th_c_m[j] as f64);
        };
        for &nj in &grid.n1[i] {
            add_neighbor(nj as usize, &mut wsum, &mut c_acc, &mut thc_acc);
        }
        for &nj in &grid.n2[i] {
            add_neighbor(nj as usize, &mut wsum, &mut c_acc, &mut thc_acc);
        }
        c_new[i] = (c_acc / wsum).clamp(0.0, 1.0) as f32;
        let advected = (thc_acc / wsum) as f32;
        let raw = advected - th_c_m[i];
        if raw.abs() > 1000.0 {
            warn_count += 1;
            if raw < warn_min {
                warn_min = raw;
            }
            if raw > warn_max {
                warn_max = raw;
            }
            warn_sum += raw as f64;
        }
        // No hard cap here; enforce global bounds after write
        thc_new[i] = th_c_m[i] + raw;
    }
    c.copy_from_slice(&c_new);
    th_c_m.copy_from_slice(&thc_new);
    // Enforce continental thickness bounds globally
    for t in th_c_m.iter_mut() {
        *t = (*t).clamp(25_000.0, 65_000.0);
    }
    if warn_count > 0 {
        let mean = warn_sum / (warn_count as f64);
        println!(
            "[advect] suspicious Δth_c count={} min={:.1} mean={:.1} max={:.1} m",
            warn_count, warn_min, mean, warn_max
        );
    }
}

/// Apply uplift to `depth_m` using the fields `C` (0..1) and `th_c_m` (m).
/// Negative uplift makes water shallower (land higher). This simply applies:
/// depth[i] += -(C[i] * th_c_m[i]).
pub fn apply_uplift_from_c_thc(depth_m: &mut [f32], c: &[f32], th_c_m: &[f32]) {
    let n = depth_m.len().min(c.len()).min(th_c_m.len());
    // Isostatic coupling: reference thickness and buoyancies differ under water vs air
    let pc = crate::PhysConsts::default();
    let th_ref_m = pc.th_ref_continental_m;
    let rho_m = pc.rho_m_kg_per_m3;
    let rho_c = pc.rho_c_kg_per_m3;
    let rho_w = pc.rho_w_kg_per_m3;
    let rho_a = pc.rho_air_kg_per_m3;
    for i in 0..n {
        let th_anom = th_c_m[i] - th_ref_m; // m; positive if thicker than reference
        let elev = -depth_m[i]; // elevation (m); >0 land, <=0 ocean
        let rho_top = if elev > 0.0 { rho_a } else { rho_w };
        let buoyancy = (rho_m - rho_c) / (rho_m - rho_top);
        let dz = -(c[i].clamp(0.0, 1.0) * buoyancy * th_anom);
        depth_m[i] = (depth_m[i] + dz).clamp(-8000.0, 8000.0);
    }
}

/// Rigid advection of continental fields under plate rotation.
/// Back-rotate each cell center by the plate's Δθ and sample nearest previous cell.
pub fn rigid_advect_c_thc(
    grid: &crate::grid::Grid,
    plates: &crate::plates::Plates,
    dt_years: f64,
    c: &mut [f32],
    th_c_m: &mut [f32],
    c_stage: &mut [f32],
    th_stage: &mut [f32],
) {
    let n = grid.cells.min(c.len()).min(th_c_m.len());
    let c_prev = c.to_owned();
    let th_prev = th_c_m.to_owned();
    let threads = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).clamp(1, 16);
    let chunk = ((n + threads - 1) / threads).max(1);
    let mut c_out: Vec<f32> = vec![0.0; n];
    let mut th_out: Vec<f32> = vec![0.0; n];
    let c_prev_ref = &c_prev;
    let th_prev_ref = &th_prev;
    std::thread::scope(|scope| {
        for t in 0..threads {
            let start = t * chunk;
            if start >= n {
                break;
            }
            let end = (start + chunk).min(n);
            let c_prev_ref = c_prev_ref;
            let th_prev_ref = th_prev_ref;
            let handle = scope.spawn(move || {
                let mut lc: Vec<f32> = vec![0.0; end - start];
                let mut lt: Vec<f32> = vec![0.0; end - start];
                for k in 0..(end - start) {
                    let i = start + k;
                    let pid = plates.plate_id.get(i).copied().unwrap_or(0) as usize;
                    let axis = [
                        plates.pole_axis[pid][0] as f64,
                        plates.pole_axis[pid][1] as f64,
                        plates.pole_axis[pid][2] as f64,
                    ];
                    let theta = -(plates.omega_rad_yr[pid] as f64) * dt_years;
                    let r = [
                        grid.pos_xyz[i][0] as f64,
                        grid.pos_xyz[i][1] as f64,
                        grid.pos_xyz[i][2] as f64,
                    ];
                    let src = crate::plates::rotate_about_axis_f64(r, axis, theta);
                    let j = nearest_from_start(grid, i, src);
                    lc[k] = c_prev_ref[j];
                    lt[k] = th_prev_ref[j];
                }
                (start, lc, lt)
            });
            let (s0, lc, lt) = handle.join().unwrap();
            c_out[s0..s0 + lc.len()].copy_from_slice(&lc);
            th_out[s0..s0 + lt.len()].copy_from_slice(&lt);
        }
    });
    c_stage[..n].copy_from_slice(&c_out[..n]);
    th_stage[..n].copy_from_slice(&th_out[..n]);
    c[..n].copy_from_slice(&c_out[..n]);
    th_c_m[..n].copy_from_slice(&th_out[..n]);
}

/// One-shot rigid remap from snapshot slices using a single composed rotation per plate.
/// Uses weighted majority (1/d) for C and masked weighted average for th_c_m over 1-ring of donor.
pub fn rigid_advect_c_thc_from(
    grid: &crate::grid::Grid,
    plates: &crate::plates::Plates,
    dt_years: f64,
    c_src: &[f32],
    th_src: &[f32],
    c_out: &mut [f32],
    th_out: &mut [f32],
) {
    let n = grid.cells.min(c_src.len()).min(th_src.len()).min(c_out.len()).min(th_out.len());
    let pc = crate::PhysConsts::default();
    let r_earth = pc.r_earth_m;
    // 1) Build forward push mapping per plate: donor j -> best target i (same plate only)
    let mut plate_to_members: std::collections::HashMap<u32, Vec<usize>> =
        std::collections::HashMap::new();
    for (i, &pid) in plates.plate_id.iter().enumerate().take(n) {
        plate_to_members.entry(u32::from(pid)).or_default().push(i);
    }
    let plate_of =
        |idx: usize| -> u32 { u32::from(plates.plate_id.get(idx).copied().unwrap_or(0)) };

    // Precompute forward rotated target positions for push mapping
    let mut tgt_pos_all: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    for (i, tp) in tgt_pos_all.iter_mut().enumerate().take(n) {
        let pid = plate_of(i) as usize;
        let axis = [
            plates.pole_axis[pid][0] as f64,
            plates.pole_axis[pid][1] as f64,
            plates.pole_axis[pid][2] as f64,
        ];
        let theta = (plates.omega_rad_yr[pid] as f64) * dt_years; // forward (push)
        let r = [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
        *tp = crate::plates::rotate_about_axis_f64(r, axis, theta);
    }

    // Build boundary rings per plate (cells within 1-ring of any cross-plate edge)
    let mut is_boundary_ring: Vec<bool> = vec![false; n];
    for (i, flag) in is_boundary_ring.iter_mut().enumerate().take(n) {
        let pid = plate_of(i);
        for &u in &grid.n1[i] {
            let j = u as usize;
            if plate_of(j) != pid {
                *flag = true;
                break;
            }
        }
    }

    let mut donor_of_target: Vec<usize> = vec![usize::MAX; n];
    let mut taken_target: Vec<bool> = vec![false; n];
    let mut ridge_leak_pre = 0usize; // pre-resolve diagnostic
    let mut holes_total = 0usize;
    let mut filled_interior = 0usize;
    let mut filled_ring = 0usize;
    for (plate_id, members) in plate_to_members.iter() {
        // Candidate pairs for this plate
        let mut pairs: Vec<(usize, usize, f64)> = Vec::with_capacity(members.len()); // (j -> i, dist)
        for &j in members {
            // Find nearest target on same plate starting at j
            let start = j;
            let i0 = nearest_from_start(grid, start, tgt_pos_all[j]);
            if plate_of(i0) != *plate_id {
                ridge_leak_pre += 1;
                continue;
            }
            let p = [
                grid.pos_xyz[i0][0] as f64,
                grid.pos_xyz[i0][1] as f64,
                grid.pos_xyz[i0][2] as f64,
            ];
            let mut d = crate::geo::great_circle_arc_len_m(tgt_pos_all[j], p, r_earth);
            // Penalize donors in boundary ring to encourage interior use first
            if is_boundary_ring[j] {
                d += 0.25 * grid.mean_cell_angle_rad() * r_earth;
            }
            pairs.push((j, i0, d));
        }
        // Sort nearest-first and assign targets uniquely
        pairs.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
        // First pass: interior donors only
        for (j, i0, _d) in pairs.iter().copied() {
            if is_boundary_ring[j] {
                continue;
            }
            if !taken_target[i0] {
                donor_of_target[i0] = j;
                taken_target[i0] = true;
            }
        }
        // Second pass: allow ring donors of same plate
        for (j, i0, _d) in pairs {
            if !is_boundary_ring[j] {
                continue;
            }
            if !taken_target[i0] {
                donor_of_target[i0] = j;
                taken_target[i0] = true;
            }
        }
        // Fill holes with nearest unused donors on same plate
        let mut used_map: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &d in donor_of_target.iter() {
            if d != usize::MAX {
                used_map.insert(d);
            }
        }
        // Prefer interior donors for hole filling
        let mut unused_interior: Vec<usize> = members
            .iter()
            .copied()
            .filter(|j| !used_map.contains(j) && !is_boundary_ring[*j])
            .collect();
        let mut unused_ring: Vec<usize> = members
            .iter()
            .copied()
            .filter(|j| !used_map.contains(j) && is_boundary_ring[*j])
            .collect();
        let hole_ids: Vec<usize> = members.iter().copied().filter(|&i| !taken_target[i]).collect();
        holes_total += hole_ids.len();
        // Greedy fill: assign each hole the nearest unused donor
        // Pass A: interior donors
        for i in hole_ids.iter().copied() {
            if taken_target[i] {
                continue;
            }
            if let Some((best_k, _best_j, _best_d)) = unused_interior
                .iter()
                .enumerate()
                .map(|(k, &j)| {
                    let p = [
                        grid.pos_xyz[i][0] as f64,
                        grid.pos_xyz[i][1] as f64,
                        grid.pos_xyz[i][2] as f64,
                    ];
                    let d = crate::geo::great_circle_arc_len_m(tgt_pos_all[j], p, r_earth);
                    (k, j, d)
                })
                .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
            {
                let j = unused_interior.remove(best_k);
                donor_of_target[i] = j;
                taken_target[i] = true;
                filled_interior += 1;
            }
        }
        // Pass B: ring donors (same plate)
        for i in hole_ids {
            if taken_target[i] {
                continue;
            }
            if let Some((best_k, _best_j, _best_d)) = unused_ring
                .iter()
                .enumerate()
                .map(|(k, &j)| {
                    let p = [
                        grid.pos_xyz[i][0] as f64,
                        grid.pos_xyz[i][1] as f64,
                        grid.pos_xyz[i][2] as f64,
                    ];
                    let d = crate::geo::great_circle_arc_len_m(tgt_pos_all[j], p, r_earth);
                    (k, j, d)
                })
                .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
            {
                let j = unused_ring.remove(best_k);
                donor_of_target[i] = j;
                taken_target[i] = true;
                filled_ring += 1;
            }
        }
    }

    // 2) Write C from bijection (no vote). Compute confidence for area correction.
    let mut conf: Vec<f32> = vec![0.0; n];
    for i in 0..n {
        let j = donor_of_target[i];
        if j == usize::MAX {
            c_out[i] = 0.0;
            th_out[i] = 0.0;
            continue;
        }
        c_out[i] = if c_src[j] > 0.5 { 1.0 } else { 0.0 };
        // Confidence: 1/d to donor + neighbors of donor on same plate
        let p_i = [grid.pos_xyz[i][0] as f64, grid.pos_xyz[i][1] as f64, grid.pos_xyz[i][2] as f64];
        let mut w_sum = 0.0;
        let mut land_w = 0.0;
        let mut th_w = 0.0;
        let mut th_acc = 0.0;
        let mut add = |jj: usize| {
            if jj >= n {
                return;
            }
            if plate_of(jj) != plate_of(j) {
                return;
            }
            let pj = [
                grid.pos_xyz[jj][0] as f64,
                grid.pos_xyz[jj][1] as f64,
                grid.pos_xyz[jj][2] as f64,
            ];
            let d = crate::geo::great_circle_arc_len_m(p_i, pj, r_earth).max(1.0);
            let w = 1.0 / d;
            w_sum += w;
            let cj = c_src[jj] > 0.5;
            if cj {
                land_w += w;
                th_acc += w * (th_src[jj] as f64);
                th_w += w;
            }
        };
        add(j);
        for &u in &grid.n1[j] {
            add(u as usize);
        }
        for &u in &grid.n2[j] {
            add(u as usize);
        }
        conf[i] = if w_sum > 0.0 { (land_w / w_sum) as f32 } else { 0.0 };
        if c_out[i] == 1.0 {
            th_out[i] = if th_w > 0.0 { (th_acc / th_w) as f32 } else { th_src[j] };
        } else {
            th_out[i] = 0.0;
        }

        // Coastal stickiness: if near previous coast and low confidence, bias to previous target label
        let c_prev_tgt = c_src[i] > 0.5;
        if (c_out[i] > 0.5) != c_prev_tgt {
            let mut coast_n1 = 0usize;
            let mut land_n1 = 0usize;
            for &u in &grid.n1[i] {
                let jj = u as usize;
                if c_src[jj] > 0.5 {
                    land_n1 += 1;
                }
                coast_n1 += 1;
            }
            let on_prev_coast = land_n1 > 0 && land_n1 < coast_n1;
            if on_prev_coast && (conf[i] as f64) < 0.65 {
                // revert to previous target label
                c_out[i] = if c_prev_tgt { 1.0 } else { 0.0 };
                if c_out[i] == 1.0 {
                    th_out[i] = th_src[i];
                } else {
                    th_out[i] = 0.0;
                }
            }
        }
    }

    // 3) Per-plate binary area correction using confidence
    let area = &grid.area;
    let mut area_before_by_plate: std::collections::HashMap<u32, f64> =
        std::collections::HashMap::new();
    let mut area_after_by_plate: std::collections::HashMap<u32, f64> =
        std::collections::HashMap::new();
    for i in 0..n {
        let pid = plate_of(i);
        *area_before_by_plate.entry(pid).or_insert(0.0) +=
            (c_src[i].clamp(0.0, 1.0) as f64) * (area[i] as f64);
        *area_after_by_plate.entry(pid).or_insert(0.0) += (c_out[i] as f64) * (area[i] as f64);
    }
    for (pid, _) in area_before_by_plate.clone() {
        let a_b = *area_before_by_plate.get(&pid).unwrap_or(&0.0);
        let a_a = *area_after_by_plate.get(&pid).unwrap_or(&0.0);
        let mut need = a_b - a_a; // +ve means need more land
        if need.abs() <= 0.0 {
            continue;
        }
        // Build candidate indices for flips sorted by confidence asc (least certain first)
        let mut cand: Vec<(usize, f64)> = Vec::new();
        for i in 0..n {
            if plate_of(i) != pid {
                continue;
            }
            let ai = area[i] as f64;
            if need > 0.0 && c_out[i] == 0.0 {
                cand.push((i, (conf[i] as f64) / ai)); // low confidence first
            } else if need < 0.0 && c_out[i] == 1.0 {
                // encode sign by negating score so sort ascending prefers low-confidence 1's
                cand.push((i, -((conf[i] as f64) / ai)));
            }
        }
        cand.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (i, score) in cand {
            if need.abs() <= 0.0 {
                break;
            }
            let ai = area[i] as f64;
            if need > 0.0 && c_out[i] == 0.0 && score >= 0.0 {
                c_out[i] = 1.0;
                need -= ai;
            } else if need < 0.0 && c_out[i] == 1.0 && score < 0.0 {
                c_out[i] = 0.0;
                th_out[i] = 0.0;
                need += ai;
            }
        }
    }

    // 4) Tiny topology cleaner then re-apply area correction
    let c_before_clean = c_out.to_vec();
    let mut c_clean = c_out.to_vec();
    for i in 0..n {
        let mut land_n1 = 0usize;
        let mut neigh = 0usize;
        for &u in &grid.n1[i] {
            let j = u as usize;
            if c_before_clean[j] > 0.5 {
                land_n1 += 1;
            }
            neigh += 1;
        }
        if c_before_clean[i] > 0.5 && land_n1 < 2 {
            c_clean[i] = 0.0;
        }
        if c_before_clean[i] <= 0.5 && land_n1 >= 4 {
            c_clean[i] = 1.0;
        }
        let _ = neigh; // silent use; thresholds are fixed
    }
    // Adjust thickness for flips introduced by cleaner
    for i in 0..n {
        if c_clean[i] > 0.5 && c_out[i] <= 0.5 {
            // became land: average neighbor thickness where land after clean or fallback to th_src[i]
            let mut acc = 0.0f64;
            let mut w = 0.0f64;
            for &u in &grid.n1[i] {
                let j = u as usize;
                if c_clean[j] > 0.5 {
                    acc += th_out[j] as f64;
                    w += 1.0;
                }
            }
            th_out[i] = if w > 0.0 { (acc / w) as f32 } else { th_src[i] };
        } else if c_clean[i] <= 0.5 && c_out[i] > 0.5 {
            th_out[i] = 0.0;
        }
    }
    c_out.copy_from_slice(&c_clean);

    // Recompute per-plate areas and correct again lightly
    area_after_by_plate.clear();
    for i in 0..n {
        let pid = plate_of(i);
        *area_after_by_plate.entry(pid).or_insert(0.0) += (c_out[i] as f64) * (area[i] as f64);
    }
    for (pid, _) in area_before_by_plate.clone() {
        let a_b = *area_before_by_plate.get(&pid).unwrap_or(&0.0);
        let a_a = *area_after_by_plate.get(&pid).unwrap_or(&0.0);
        let mut need = a_b - a_a;
        if need.abs() <= 0.0 {
            continue;
        }
        let mut cand: Vec<(usize, f64)> = Vec::new();
        for i in 0..n {
            if plate_of(i) != pid {
                continue;
            }
            let ai = area[i] as f64;
            if need > 0.0 && c_out[i] == 0.0 {
                cand.push((i, (conf[i] as f64) / ai));
            } else if need < 0.0 && c_out[i] == 1.0 {
                cand.push((i, -((conf[i] as f64) / ai)));
            }
        }
        cand.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (i, score) in cand {
            if need.abs() <= 0.0 {
                break;
            }
            let ai = area[i] as f64;
            if need > 0.0 && c_out[i] == 0.0 && score >= 0.0 {
                c_out[i] = 1.0;
                need -= ai;
            } else if need < 0.0 && c_out[i] == 1.0 && score < 0.0 {
                c_out[i] = 0.0;
                th_out[i] = 0.0;
                need += ai;
            }
        }
    }

    // 4) Report diagnostics
    // Final diagnostics: collisions, ridge_leak after assignment, and hole fill breakdown
    let mut shared2 = 0u32; // Should be zero with strict bijection
    let mut used_donors = std::collections::HashSet::new();
    for &d in donor_of_target.iter().take(n) {
        if d != usize::MAX && !used_donors.insert(d) {
            shared2 += 1;
        }
    }
    let mut ridge_leak_final = 0usize;
    for (i, &d) in donor_of_target.iter().enumerate().take(n) {
        if d != usize::MAX && plate_of(d) != plate_of(i) {
            ridge_leak_final += 1;
        }
    }
    println!("[rigid][nn] collisions_final={} ridge_leak={} holes_filled={} filled_int={} filled_ring={} | leak_pre={}", shared2, ridge_leak_final, holes_total, filled_interior, filled_ring, ridge_leak_pre);
}

/// Greedy hill-climb to nearest cell center to `target` starting from `start_i`.
fn nearest_from_start(grid: &crate::grid::Grid, mut idx: usize, target: [f64; 3]) -> usize {
    let mut best = dot3_f64(
        [grid.pos_xyz[idx][0] as f64, grid.pos_xyz[idx][1] as f64, grid.pos_xyz[idx][2] as f64],
        target,
    );
    loop {
        let mut improved = false;
        for &n1 in &grid.n1[idx] {
            let j = n1 as usize;
            let d = dot3_f64(
                [grid.pos_xyz[j][0] as f64, grid.pos_xyz[j][1] as f64, grid.pos_xyz[j][2] as f64],
                target,
            );
            if d > best || ((d - best).abs() <= f64::EPSILON && j < idx) {
                best = d;
                idx = j;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    idx
}

#[inline]
fn dot3_f64(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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
