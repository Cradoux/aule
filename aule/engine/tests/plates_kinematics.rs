use engine as crate_engine;

#[test]
fn rodrigues_roundtrip_small_error() {
    let axis = [0.0f32, 0.0, 1.0];
    let r = [1.0f32, 0.0, 0.0];
    let theta = 0.25f32;
    let r1 = crate_engine::plates::rotate_point(axis, theta, r);
    let r2 = crate_engine::plates::rotate_point(axis, -theta, r1);
    let err = ((r2[0] as f64 - r[0] as f64).powi(2)
        + (r2[1] as f64 - r[1] as f64).powi(2)
        + (r2[2] as f64 - r[2] as f64).powi(2))
    .sqrt();
    assert!(err < 1e-6, "roundtrip err = {}", err);
}
