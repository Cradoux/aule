//! Small utilities for sanitizing fields and safety checks.

/// Clamp and scrub elevation array to safe display range.
pub fn sanitize_elevation(h: &mut [f32]) {
    for v in h.iter_mut() {
        if !v.is_finite() {
            *v = 0.0;
        } else {
            // Land ~[-inf, +9km], ocean ~[-11km, 0]; clamp to visible range
            *v = v.clamp(-11_000.0, 9_000.0);
        }
    }
}

/// Smoothly saturate a value toward ±cap using a tanh function (f64 variant).
///
/// Returns `cap * tanh(x / cap)`. For |x| << cap, returns ~x; for large |x|, asymptotes to ±cap.
#[inline]
pub fn soft_cap_f64(x: f64, cap: f64) -> f64 {
    if cap <= 0.0 {
        return 0.0;
    }
    let a = x / cap;
    cap * a.tanh()
}

/// Smoothly saturate a value toward ±cap using a tanh function (f32 variant).
#[inline]
pub fn soft_cap_f32(x: f32, cap: f32) -> f32 {
    if cap <= 0.0 {
        return 0.0;
    }
    let a = x / cap;
    cap * a.tanh()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn soft_cap_behaves_reasonably_f64() {
        assert_eq!(soft_cap_f64(0.0, 1.0), 0.0);
        let near = soft_cap_f64(0.1, 10.0);
        assert!((near - 0.1).abs() < 5e-6);
        let hi = soft_cap_f64(1000.0, 10.0);
        assert!(hi <= 10.0 && hi > 9.0);
        let lo = soft_cap_f64(-1000.0, 10.0);
        assert!((-10.0..-9.0).contains(&lo));
    }

    #[test]
    fn soft_cap_behaves_reasonably_f32() {
        assert_eq!(soft_cap_f32(0.0, 1.0), 0.0);
        let near = soft_cap_f32(0.1, 10.0);
        assert!((near - 0.1).abs() < 1e-5);
        let hi = soft_cap_f32(1000.0, 10.0);
        assert!(hi <= 10.0 && hi > 9.0);
        let lo = soft_cap_f32(-1000.0, 10.0);
        assert!((-10.0f32..-9.0f32).contains(&lo));
    }
}
