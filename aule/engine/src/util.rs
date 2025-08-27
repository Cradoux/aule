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
