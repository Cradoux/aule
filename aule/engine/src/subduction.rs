//! Subduction bands from convergent boundaries and bathymetry adjustments (CPU).
//!
//! Modeling notes:
//! - Seeds are picked from convergent boundary edges and split into subducting vs overriding sides
//!   by age (older plate subducts; tie-break by plate id). This is a heuristic; real polarity may
//!   depend on buoyancy (age), trench rollback, and slab geometry.
//! - Band geometries are distance-thresholded within plate domains to avoid crossing plate labels.
//! - Bathymetry edits are idempotent per call: we recompute the Parsons–Sclater baseline from age
//!   then add band deltas. Callers should ensure they apply this after baseline age→depth mapping.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::boundaries::Boundaries;
use crate::geo;
use crate::grid::Grid;

const KM: f64 = 1000.0;
const RADIUS_M: f64 = 6_371_000.0;

/// Boolean per-cell masks for each subduction-associated band.
#[derive(Default)]
pub struct SubductionMasks {
    /// Cells within trench band on the subducting plate.
    pub trench: Vec<bool>,
    /// Cells within magmatic arc band on the overriding plate.
    pub arc: Vec<bool>,
    /// Cells within back-arc band on the overriding plate.
    pub backarc: Vec<bool>,
}

/// User parameters that control band geometry and bathymetric edits.
#[derive(Clone, Copy)]
pub struct SubductionParams {
    /// Convergence threshold in m/yr (reused from classification scale).
    pub tau_conv_m_per_yr: f64,
    /// Half-width of trench band (km) from subducting seeds.
    pub trench_half_width_km: f64,
    /// Offset of peak arc band from overriding seeds (km).
    pub arc_offset_km: f64,
    /// Half-width of arc band (km).
    pub arc_half_width_km: f64,
    /// Back-arc band width (km) immediately behind arc.
    pub backarc_width_km: f64,
    /// Trench bathymetry delta (m, positive deepens).
    pub trench_deepen_m: f32,
    /// Arc bathymetry delta (m, negative uplifts/shallows).
    pub arc_uplift_m: f32,
    /// Back-arc bathymetry delta (m, negative uplifts/shallows).
    pub backarc_uplift_m: f32,
    /// Absolute oceanward offset applied to the hinge (meters). Default: 0.
    pub rollback_offset_m: f64,
    /// If > 0, can be used by viewer stepper: offset += rate * dt_myr (not applied here).
    pub rollback_rate_km_per_myr: f64,
    /// If true, back-arc is deepened instead of uplifted (extension mode).
    pub backarc_extension_mode: bool,
    /// Depth added inside back-arc band in extension mode (positive down).
    pub backarc_extension_deepen_m: f32,
    /// Continental fraction threshold to treat a cell as continental (0..1).
    pub continent_c_min: f32,
}

/// Counts for each band.
pub struct SubductionStats {
    /// Number of trench cells.
    pub trench_cells: u32,
    /// Number of arc cells.
    pub arc_cells: u32,
    /// Number of back-arc cells.
    pub backarc_cells: u32,
}

/// Result bundle: masks and summary stats.
pub struct SubductionResult {
    /// Per-band boolean masks per cell.
    pub masks: SubductionMasks,
    /// Summary statistics for band sizes.
    pub stats: SubductionStats,
}

#[derive(Copy, Clone, Debug)]
struct QItem {
    dist_m: f64,
    cell: u32,
}
impl Eq for QItem {}
impl PartialEq for QItem {
    fn eq(&self, other: &Self) -> bool {
        self.dist_m.eq(&other.dist_m) && self.cell == other.cell
    }
}
impl PartialOrd for QItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for QItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // min-heap by distance
        match self.dist_m.partial_cmp(&other.dist_m) {
            Some(Ordering::Less) => Ordering::Greater,
            Some(Ordering::Greater) => Ordering::Less,
            Some(Ordering::Equal) => self.cell.cmp(&other.cell).reverse(),
            None => Ordering::Equal,
        }
    }
}

/// Compute subduction bands and adjust bathymetry in-place (idempotent).
pub fn apply_subduction(
    grid: &Grid,
    boundaries: &Boundaries,
    plate_id: &[u16],
    age_myr: &[f32],
    v_en: &[[f32; 2]],
    depth_m: &mut [f32],
    params: SubductionParams,
    c_opt: Option<&[f32]>,
) -> SubductionResult {
    assert_eq!(plate_id.len(), grid.cells);
    assert_eq!(age_myr.len(), grid.cells);
    assert_eq!(v_en.len(), grid.cells);
    assert_eq!(depth_m.len(), grid.cells);

    // Seeds from convergent edges
    let mut sub_seeds: Vec<u32> = Vec::new();
    let mut over_seeds: Vec<u32> = Vec::new();
    for ek in &boundaries.edge_kin {
        if ek.class as u8 != 2 {
            continue;
        } // convergent
        let u = ek.u;
        let v = ek.v;
        let au = age_myr[u as usize] as f64;
        let av = age_myr[v as usize] as f64;
        let (s, o) = if (au > av) || (au == av && plate_id[u as usize] > plate_id[v as usize]) {
            (u, v)
        } else {
            (v, u)
        };
        sub_seeds.push(s);
        over_seeds.push(o);
    }

    // Precompute unit positions
    let mut pos_unit: Vec<[f64; 3]> = Vec::with_capacity(grid.cells);
    for p in &grid.pos_xyz {
        pos_unit.push(geo::normalize([p[0] as f64, p[1] as f64, p[2] as f64]));
    }

    // Multi-source Dijkstra helper restricted to plate domains
    let mut dist_sub: Vec<f64> = vec![f64::INFINITY; grid.cells];
    let mut dist_over: Vec<f64> = vec![f64::INFINITY; grid.cells];

    let mut heap: BinaryHeap<QItem> = BinaryHeap::new();

    // Subducting side
    for &s in &sub_seeds {
        dist_sub[s as usize] = 0.0;
        heap.push(QItem { dist_m: 0.0, cell: s });
    }
    while let Some(QItem { dist_m, cell }) = heap.pop() {
        let u = cell as usize;
        if dist_m > dist_sub[u] {
            continue;
        }
        let pid = plate_id[u];
        for &vn in &grid.n1[u] {
            let v = vn as usize;
            if plate_id[v] != pid {
                continue;
            }
            let s_m = geo::great_circle_arc_len_m(pos_unit[u], pos_unit[v], RADIUS_M);
            let nd = dist_m + s_m;
            if nd < dist_sub[v] {
                dist_sub[v] = nd;
                heap.push(QItem { dist_m: nd, cell: v as u32 });
            }
        }
    }
    // Overriding side
    heap.clear();
    for &s in &over_seeds {
        dist_over[s as usize] = 0.0;
        heap.push(QItem { dist_m: 0.0, cell: s });
    }
    while let Some(QItem { dist_m, cell }) = heap.pop() {
        let u = cell as usize;
        if dist_m > dist_over[u] {
            continue;
        }
        let pid = plate_id[u];
        for &vn in &grid.n1[u] {
            let v = vn as usize;
            if plate_id[v] != pid {
                continue;
            }
            let s_m = geo::great_circle_arc_len_m(pos_unit[u], pos_unit[v], RADIUS_M);
            let nd = dist_m + s_m;
            if nd < dist_over[v] {
                dist_over[v] = nd;
                heap.push(QItem { dist_m: nd, cell: v as u32 });
            }
        }
    }

    // Thresholds
    let trench_hw_m = params.trench_half_width_km * KM;
    let arc_off_m = params.arc_offset_km * KM;
    let arc_hw_m = params.arc_half_width_km * KM;
    let backarc_w_m = params.backarc_width_km * KM;

    let mut masks = SubductionMasks {
        trench: vec![false; grid.cells],
        arc: vec![false; grid.cells],
        backarc: vec![false; grid.cells],
    };

    let mut stats = SubductionStats { trench_cells: 0, arc_cells: 0, backarc_cells: 0 };

    for i in 0..grid.cells {
        let is_trench = dist_sub[i].is_finite() && dist_sub[i] <= trench_hw_m;
        let mut is_arc = false;
        let mut is_back = false;
        if dist_over[i].is_finite() {
            // Apply rollback offset on overriding distance
            let d_eff = (dist_over[i] - params.rollback_offset_m).max(0.0);
            is_arc = (d_eff - arc_off_m).abs() <= arc_hw_m;
            is_back = !is_arc
                && d_eff >= (arc_off_m + arc_hw_m)
                && d_eff <= (arc_off_m + arc_hw_m + backarc_w_m);
        }
        if is_trench {
            masks.trench[i] = true;
            stats.trench_cells += 1;
        } else if is_arc {
            masks.arc[i] = true;
            stats.arc_cells += 1;
        } else if is_back {
            masks.backarc[i] = true;
            stats.backarc_cells += 1;
        }
    }

    // Idempotent bathy edits: recompute baseline from age curve, then add deltas
    const D0: f64 = 2600.0;
    const A_COEF: f64 = 350.0;
    const B_COEF: f64 = 0.0;
    for i in 0..grid.cells {
        let mut delta: f32 = 0.0;
        if masks.trench[i] {
            // Gate trench deepening by continental fraction if provided
            let is_cont =
                c_opt.and_then(|c| c.get(i)).map(|&v| v >= params.continent_c_min).unwrap_or(false);
            if !is_cont {
                delta += params.trench_deepen_m;
            }
        }
        if masks.arc[i] {
            delta += params.arc_uplift_m;
        }
        if masks.backarc[i] {
            if params.backarc_extension_mode {
                delta += params.backarc_extension_deepen_m;
            } else {
                delta += params.backarc_uplift_m;
            }
        }
        if delta != 0.0 {
            let mut base = crate::age::depth_from_age(age_myr[i] as f64, D0, A_COEF, B_COEF) as f32;
            if !base.is_finite() {
                base = 6000.0;
            }
            base = base.clamp(0.0, 6000.0);
            depth_m[i] = base + delta;
            // NOTE: This overwrites any prior tectonic edit at `i` in the same step. If multiple
            // processes should superpose (e.g., transforms), apply them consistently after all
            // baselines or accumulate deltas before a single write.
        }
    }

    SubductionResult { masks, stats }
}
