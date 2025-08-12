use engine::flexure::{alpha, d_from_te, line_load_rhs, solve_flexure_1d, w_line_analytic};

#[test]
fn alpha_monotonicity() {
    let k = 2.5e5; // choose k so that alpha is O(60–100 km)
    let d1 = d_from_te(70e9, 0.25, 10_000.0);
    let d2 = d_from_te(70e9, 0.25, 20_000.0);
    let a1 = alpha(d1, k);
    let a2 = alpha(d2, k);
    assert!(a2 > a1);
}

#[test]
fn numeric_matches_analytic_line_load() {
    // Parameters (SI)
    let e = 70e9;
    let nu = 0.25;
    let te = 25_000.0;
    let k = 2.5e5; // larger alpha for adequate grid resolution
    let p = 1.0e12;
    let d = d_from_te(e, nu, te);
    let a = alpha(d, k);
    // Grid
    let dx = 5_000.0;
    let n = 1025usize; // ~5,120 km
    let q = line_load_rhs(n, dx, p);
    let sol = solve_flexure_1d(&q, dx, d, k, 1e-8, 20_000);
    assert!(sol.iters < 5_000);
    // Compare RMS in interior |x| <= 4α
    let i0 = n / 2;
    let xmax = 4.0 * a;
    let imax = ((xmax / dx).floor() as isize).min(i0 as isize);
    let mut sum2 = 0.0;
    let mut sum2_ref = 0.0;
    let mut count = 0.0;
    for di in -(imax as isize)..=(imax as isize) {
        let i = (i0 as isize + di) as usize;
        let x = (di as f64) * dx;
        let wa = w_line_analytic(x, p, d, k);
        let wn = sol.w[i];
        sum2 += (wn - wa) * (wn - wa);
        sum2_ref += wa * wa;
        count += 1.0;
    }
    let rms = (sum2 / count).sqrt();
    let ref_rms = (sum2_ref / count).sqrt();
    assert!(rms / ref_rms <= 0.05, "rms/peak = {}", rms / ref_rms);
}
