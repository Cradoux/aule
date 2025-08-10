use engine::grid::tile::Tiling;
use engine::grid::Grid;

#[test]
fn tiler_trivial_covers_all() {
    let g = Grid::new(16);
    let t = Tiling::new(&g, 2, 8192);
    // Trivial tiler currently produces one tile that owns all cells.
    assert_eq!(t.tiles.len(), 1);
    assert_eq!(t.owner.len(), g.cells);
    assert!(t.tiles[0].halo.is_empty());
    assert_eq!(t.tiles[0].interior.len(), g.cells);
    // owner rule
    for &own in &t.owner {
        assert_eq!(own, 0);
    }
}
