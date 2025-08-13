use viewer::colormap::{parse_pal, sample_linear_srgb, HYPS_DEFAULT_STR};

#[test]
fn parse_pal_increasing_and_minmax() {
    let pal = parse_pal(HYPS_DEFAULT_STR).expect("parse pal");
    assert!(pal.stops.len() >= 2);
    for i in 1..pal.stops.len() {
        assert!(
            pal.stops[i].v.partial_cmp(&pal.stops[i - 1].v) == Some(std::cmp::Ordering::Greater)
        );
    }
    assert!(pal.vmax > pal.vmin);
}

#[test]
fn interpolation_midpoint_linear_rgb() {
    // Simple two-stop palette: black at 0, white at 100
    let src = "0 #000000\n100 #FFFFFF\n";
    let pal = parse_pal(src).unwrap();
    let rgb = sample_linear_srgb(&pal, 50.0);
    // Expect ~#BBBBBB (≈187) in sRGB due to gamma (not 128); accept ±1
    assert!((rgb[0] as i32 - 187).abs() <= 1);
    assert!((rgb[1] as i32 - 187).abs() <= 1);
    assert!((rgb[2] as i32 - 187).abs() <= 1);
}

#[test]
fn clamping_at_ends() {
    let src = "0 #112233\n10 #445566\n";
    let pal = parse_pal(src).unwrap();
    let a = sample_linear_srgb(&pal, -5.0);
    let b = sample_linear_srgb(&pal, 15.0);
    assert_eq!(a, [0x11, 0x22, 0x33]);
    assert_eq!(b, [0x44, 0x55, 0x66]);
}
