use engine::hypsometry::histogram_area_weighted;

#[test]
fn area_conservation_and_land_fraction() {
    let depth = vec![100.0, -50.0, 0.0, 200.0, -300.0];
    let area = vec![1.0e6, 2.0e6, 1.5e6, 3.0e6, 2.5e6];
    let bins = 10usize;
    let (centers, areas) = histogram_area_weighted(&depth, &area, bins, -400.0, 300.0);
    assert_eq!(centers.len(), bins);
    assert_eq!(areas.len(), bins);
    let total_area_hist: f64 = areas.iter().sum();
    let total_area_direct: f64 = area.iter().map(|&a| a as f64).sum();
    let rel = ((total_area_hist - total_area_direct) / total_area_direct).abs();
    assert!(rel < 1e-12, "area rel err {}", rel);

    let land_direct: f64 =
        depth.iter().zip(area.iter()).filter(|(&d, _)| d <= 0.0).map(|(_, &a)| a as f64).sum();
    let land_hist: f64 =
        areas.iter().zip(centers.iter()).filter(|(_, &c)| c <= 0.0).map(|(a, _)| *a).sum();
    let err = (land_hist - land_direct).abs() / land_direct.max(1.0);
    assert!(err < 1e-12, "land frac err {}", err);

    // determinism
    let (c2, a2) = histogram_area_weighted(&depth, &area, bins, -400.0, 300.0);
    assert_eq!(centers, c2);
    assert_eq!(areas, a2);
}
