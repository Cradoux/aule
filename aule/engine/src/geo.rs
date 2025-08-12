#![allow(clippy::many_single_char_names)]

/// Dot product of 3D vectors.
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of 3D vectors.
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Euclidean norm of a 3D vector.
#[inline]
pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Normalize a 3D vector (returns zero if input is zero).
#[inline]
pub fn normalize(mut a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n > 0.0 {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
    a
}

/// Great-circle arc length (meters) between two unit vectors on a sphere of radius `r_m`.
#[inline]
pub fn great_circle_arc_len_m(a_unit: [f64; 3], b_unit: [f64; 3], r_m: f64) -> f64 {
    // angle = atan2(|a×b|, a·b)
    let c = cross(a_unit, b_unit);
    let s = norm(c);
    let d = dot(a_unit, b_unit).clamp(-1.0, 1.0);
    r_m * s.atan2(d).abs()
}

/// Local tangent basis at unit position r_hat: (east, north).
#[inline]
pub fn local_basis(r_hat: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Use Z as reference; if near-parallel, fall back to Y
    let z = [0.0_f64, 0.0, 1.0];
    let mut e = cross(z, r_hat);
    if norm(e) < 1e-12 {
        let y = [0.0_f64, 1.0, 0.0];
        e = cross(y, r_hat);
    }
    let east = normalize(e);
    let north = normalize(cross(r_hat, east));
    (east, north)
}

/// Convert local (east,north) components at r̂ into a world-space vector (m/yr).
#[inline]
pub fn en_to_world(r_hat: [f64; 3], v_en_m_per_yr: [f32; 2]) -> [f64; 3] {
    let (e, n) = local_basis(r_hat);
    let ve = v_en_m_per_yr[0] as f64;
    let vn = v_en_m_per_yr[1] as f64;
    [e[0] * ve + n[0] * vn, e[1] * ve + n[1] * vn, e[2] * ve + n[2] * vn]
}
