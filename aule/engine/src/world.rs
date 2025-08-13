//! World state container and constructors (T-400).

use crate::{boundaries::Boundaries, grid::Grid, plates::Plates};

/// Simulation clock information.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Clock {
    /// Current simulation time in Myr.
    pub t_myr: f64,
    /// Step index (starts at 0, increments per step).
    pub step_idx: u64,
}

/// The complete world state required for stepping.
pub struct World {
    /// Geodesic grid definition.
    pub grid: Grid,
    /// Plate assignment and Euler poles.
    pub plates: Plates,
    /// Classified boundaries for current state.
    pub boundaries: Boundaries,
    /// Age field in Myr, length = grid.cells.
    pub age_myr: Vec<f32>,
    /// Bathymetry depth in meters (+ down), length = grid.cells.
    pub depth_m: Vec<f32>,
    /// Per-cell velocities (east,north) in m/yr.
    pub v_en: Vec<[f32; 2]>,
    /// Simulation clock.
    pub clock: Clock,
    /// Reference sea-level volume and area captured after baseline age→depth
    pub sea_level_ref: Option<SeaLevelRef>,
    /// Precomputed per-cell areas in m^2 on the sphere (4πR^2 scaled)
    pub area_m2: Vec<f32>,
}

/// Reference sea-level bookkeeping
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeaLevelRef {
    /// Target ocean volume at reference state (m^3)
    pub volume_m3: f64,
    /// Ocean area used at reference state (m^2)
    pub ocean_area_m2: f64,
}

impl World {
    /// Construct a world with deterministic plates and zeroed ages/depths.
    pub fn new(f: u32, num_plates: u32, seed: u64) -> Self {
        let grid = Grid::new(f);
        let plates = Plates::new(&grid, num_plates, seed);
        let v_en = plates.vel_en.clone();
        let boundaries =
            crate::boundaries::Boundaries::classify(&grid, &plates.plate_id, &v_en, 0.005);
        let age_myr = vec![0.0f32; grid.cells];
        let depth_m = vec![0.0f32; grid.cells];
        let clock = Clock { t_myr: 0.0, step_idx: 0 };
        // Precompute area in m^2 using Earth radius
        const R_EARTH_M: f64 = 6_371_000.0;
        let scale = 4.0 * std::f64::consts::PI * R_EARTH_M * R_EARTH_M;
        let mut area_m2: Vec<f32> = Vec::with_capacity(grid.cells);
        for &a in &grid.area {
            area_m2.push((a as f64 * scale) as f32);
        }
        Self {
            grid,
            plates,
            boundaries,
            age_myr,
            depth_m,
            v_en,
            clock,
            sea_level_ref: None,
            area_m2,
        }
    }
}
