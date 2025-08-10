use engine::grid::{cache::GridError, Grid};

#[test]
fn grid_constructs() {
    let g = Grid::new(0);
    assert_eq!(g.cells, 12);
    assert_eq!(g.pos_xyz.len(), g.cells);
    assert_eq!(g.latlon.len(), g.cells);
    assert_eq!(g.area.len(), g.cells);
    assert_eq!(g.n1.len(), g.cells);
    assert_eq!(g.n2.len(), g.cells);
}

#[test]
fn cache_round_trip_minimal() -> Result<(), GridError> {
    let g = Grid::new(0);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("grid.cache");
    g.save_cache(&path)?;
    let g2 = Grid::load_cache(&path)?;
    assert_eq!(g.cells, g2.cells);
    assert_eq!(g.pos_xyz, g2.pos_xyz);
    assert_eq!(g.latlon, g2.latlon);
    assert_eq!(g.area, g2.area);
    assert_eq!(g.n1, g2.n1);
    assert_eq!(g.n2, g2.n2);
    Ok(())
}
