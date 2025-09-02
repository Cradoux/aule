//! 1D elastic plate flexure on a Winkler foundation (CPU reference).
//!
//! Equation: D w''''(x) + k w(x) = q(x)
//! - D: plate rigidity (Pa·m^3)
//! - k: foundation stiffness (N/m^3)
//! - q: load per unit length (N/m)
//!
//! Analytic line load solution for q(x) = P δ(x):
//! w(x) = (P / (2 k α)) exp(-|x|/α) [cos(|x|/α) + sin(|x|/α)]

/// Compute plate rigidity from elastic thickness.
pub fn d_from_te(e_pa: f64, nu: f64, te_m: f64) -> f64 {
    let denom = 12.0 * (1.0 - nu * nu);
    e_pa * te_m.powi(3) / denom
}

/// Flexural parameter α = (4D/k)^(1/4).
pub fn alpha(d: f64, k: f64) -> f64 {
    ((4.0 * d) / k).powf(0.25)
}

/// Analytic deflection for a line load on an infinite plate.
pub fn w_line_analytic(x: f64, p_n_per_m: f64, d: f64, k: f64) -> f64 {
    let a = alpha(d, k);
    let s = (x.abs()) / a;
    let pref = p_n_per_m / (2.0 * k * a);
    pref * (-s).exp() * (s.cos() + s.sin())
}

/// Compute a simple Te field (m) from age (Myr) and continental fraction C (0..1).
///
/// Policy:
/// - Oceanic Te grows with age, from 5 km at 0 Myr to 40 km at ≥120 Myr (logistic cap).
/// - Continental Te uses a higher baseline: 20–70 km depending on th_c and C.
/// - Blend by C.
pub fn compute_te_field(age_myr: &[f32], c: &[f32], th_c_m: &[f32], out_te_m: &mut [f32]) {
    let n = out_te_m.len().min(age_myr.len()).min(c.len()).min(th_c_m.len());
    for i in 0..n {
        let age = age_myr[i] as f64;
        let ci = c[i].clamp(0.0, 1.0) as f64;
        // Oceanic: logistic growth with age
        let te_ocean = {
            let te_min = 5_000.0f64;
            let te_max = 40_000.0f64;
            let t50 = 30.0f64; // age at half-rise (Myr)
            let k = 0.08f64; // steepness
            te_min + (te_max - te_min) / (1.0 + (-k * (age - t50)).exp())
        };
        // Continental: scale with thickness (20–70 km across 20–70 km crust)
        let thc = th_c_m[i] as f64;
        let te_cont = 20_000.0 + (thc - 20_000.0).clamp(0.0, 50_000.0);
        let te = (1.0 - ci) * te_ocean + ci * te_cont;
        out_te_m[i] = te as f32;
    }
}

/// Result of the CG solve for 1D flexure.
#[derive(Debug, Clone)]
pub struct FlexSolve {
    /// Deflection (m), length N.
    pub w: Vec<f64>,
    /// Iterations used.
    pub iters: usize,
    /// Final residual 2-norm.
    pub resid: f64,
}

/// Apply bi-harmonic 5-point stencil to vector `v` into `out` (zero BC beyond ends).
fn apply_biharmonic_5(v: &[f64], out: &mut [f64], inv_dx4: f64) {
    let n = v.len();
    out.fill(0.0);
    if n == 0 {
        return;
    }
    for i in 0..n {
        let w_im2 = if i >= 2 { v[i - 2] } else { 0.0 };
        let w_im1 = if i >= 1 { v[i - 1] } else { 0.0 };
        let w_ip1 = if i + 1 < n { v[i + 1] } else { 0.0 };
        let w_ip2 = if i + 2 < n { v[i + 2] } else { 0.0 };
        let val = (w_im2 - 4.0 * w_im1 + 6.0 * v[i] - 4.0 * w_ip1 + w_ip2) * inv_dx4;
        out[i] = val;
    }
}

/// Conjugate Gradient solve of (D*bi4 + k*I) w = q, with homogeneous Dirichlet at ends via zero stencil.
pub fn solve_flexure_1d(
    q: &[f64],
    dx: f64,
    d: f64,
    k: f64,
    tol: f64,
    max_iter: usize,
) -> FlexSolve {
    let n = q.len();
    let inv_dx4 = 1.0 / (dx * dx * dx * dx);
    let mut w = vec![0.0f64; n];
    let mut r = q.to_vec();
    // r = q - A w (w=0 initially)
    let mut p = r.clone();
    let mut ap = vec![0.0f64; n];
    let mut tmp = vec![0.0f64; n];
    let mut rsold = r.iter().map(|x| x * x).sum::<f64>();
    let mut resid = rsold.sqrt();
    let mut iters = 0usize;
    if resid <= tol {
        return FlexSolve { w, iters, resid };
    }
    for it in 0..max_iter {
        // ap = D*bi4(p) + k*p
        apply_biharmonic_5(&p, &mut tmp, inv_dx4);
        for i in 0..n {
            ap[i] = d * tmp[i] + k * p[i];
        }
        let denom = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum::<f64>();
        if denom == 0.0 {
            break;
        }
        let alpha_cg = rsold / denom;
        for i in 0..n {
            w[i] += alpha_cg * p[i];
        }
        for i in 0..n {
            r[i] -= alpha_cg * ap[i];
        }
        let rsnew = r.iter().map(|x| x * x).sum::<f64>();
        resid = rsnew.sqrt();
        iters = it + 1;
        if resid <= tol {
            break;
        }
        let beta = rsnew / rsold;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    FlexSolve { w, iters, resid }
}

/// Build a symmetric grid RHS for a single line load centered at index i0.
pub fn line_load_rhs(n: usize, dx: f64, p_n_per_m: f64) -> Vec<f64> {
    let mut q = vec![0.0f64; n];
    if n == 0 {
        return q;
    }
    let i0 = n / 2;
    q[i0] = p_n_per_m / dx;
    q
}
