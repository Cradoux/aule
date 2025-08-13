//! T-455: tiling plan tests (CPU-only staging).

use engine::{grid::Grid, tileplan};

fn build_plan_default(grid: &Grid) -> tileplan::TilePlan {
    let p = tileplan::TileBuildParams {
        target_cells_per_tile: 4096,
        halo_rings: 1,
        max_tile_cells: 32_768,
    };
    tileplan::build_tileplan(grid, p).expect("tileplan build")
}

#[test]
fn halo_completeness_one_ring() {
    let g = Grid::new(64);
    let plan = build_plan_default(&g);
    tileplan::validate_tileplan(&g, &plan, 1).expect("valid halo");
}

#[test]
fn determinism_same_input_same_plan() {
    let g = Grid::new(64);
    let p0 = build_plan_default(&g);
    let p1 = build_plan_default(&g);
    assert_eq!(p0.halo_rings, p1.halo_rings);
    assert_eq!(p0.max_tile_cells, p1.max_tile_cells);
    assert_eq!(p0.tiles.len(), p1.tiles.len());
    for (a, b) in p0.tiles.iter().zip(p1.tiles.iter()) {
        assert_eq!(a.interior, b.interior);
        assert_eq!(a.halo, b.halo);
    }
}

#[test]
fn memory_guard_too_small() {
    let g = Grid::new(32);
    let p =
        tileplan::TileBuildParams { target_cells_per_tile: 64, halo_rings: 1, max_tile_cells: 10 };
    let res = tileplan::build_tileplan(&g, p);
    assert!(matches!(res, Err(tileplan::TileError::TooSmall(_, 10))));
}
