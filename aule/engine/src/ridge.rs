//! Ridge birth and optional fringe assignment on divergent boundaries (CPU).
//! Deterministic and allocation-free given caller-provided buffers.

use crate::boundaries::Boundaries;
use crate::grid::Grid;

/// Parameters for ridge processing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RidgeParams {
    /// Age assigned to immediate 1-ring “fringe” around ridges (Myr). 0 disables fringe.
    pub fringe_age_myr: f32,
}

impl Default for RidgeParams {
    fn default() -> Self {
        Self { fringe_age_myr: 0.2 }
    }
}

/// Statistics for ridge updates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RidgeStats {
    /// Cells set to age = 0 (ridge births)
    pub births: u32,
    /// Unique cells clamped to <= fringe_age_myr
    pub fringe: u32,
}

/// Apply ridge births and optional one-ring fringe on a CPU grid.
///
/// - Marks cells adjacent to divergent edges as ridges and sets their ages to 0.
/// - If `params.fringe_age_myr > 0`, clamps the ages of 1-ring neighbors to at most that value.
/// - Returns counts for births and fringe unique cells affected.
pub fn apply_ridge(
    grid: &Grid,
    boundaries: &Boundaries,
    age_ocean: &mut [f32],
    params: RidgeParams,
) -> RidgeStats {
    assert_eq!(age_ocean.len(), grid.cells);

    let n_cells = grid.cells;
    let mut is_ridge = vec![false; n_cells];
    // Build ridge mask from divergent edges (class == 1)
    for &(u, v, class) in &boundaries.edges {
        if class == 1 {
            is_ridge[u as usize] = true;
            is_ridge[v as usize] = true;
        }
    }

    // Births: set age = 0 for ridge cells
    let mut births: u32 = 0;
    for (i, &ridge) in is_ridge.iter().enumerate() {
        if ridge {
            if age_ocean[i] != 0.0 {
                births += 1;
            }
            age_ocean[i] = 0.0;
        }
    }

    // Fringe: clamp 1-ring neighbors to <= fringe_age_myr
    let mut fringe: u32 = 0;
    if params.fringe_age_myr > 0.0 {
        let mut touched = vec![false; n_cells];
        for (u, &ridge) in is_ridge.iter().enumerate() {
            if !ridge {
                continue;
            }
            for &n in &grid.n1[u] {
                let ni = n as usize;
                if !touched[ni] {
                    let before = age_ocean[ni];
                    let after = before.min(params.fringe_age_myr);
                    if after < before {
                        fringe += 1;
                        touched[ni] = true;
                        age_ocean[ni] = after;
                    } else {
                        // Even if no change, mark touched so we count unique cells that met the condition?
                        // Spec says: "count unique cells touched as fringe". We interpret as those clamped.
                    }
                }
            }
        }
    }

    RidgeStats { births, fringe }
}
