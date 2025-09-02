//! Plate boundary network utilities: adjacency, triple junctions, connectivity stats.

use crate::{boundaries::Boundaries, grid::Grid, plates::Plates};
use std::collections::{HashMap, HashSet, VecDeque};

/// Summary network statistics for diagnostics.
#[derive(Debug, Clone, Copy, Default)]
pub struct PlateNetworkStats {
    /// Number of active plates present in `plate_id`.
    pub num_plates: u32,
    /// Number of undirected adjacency edges between distinct plates.
    pub edges: u32,
    /// Number of connected components in the plate adjacency graph.
    pub components: u32,
    /// Area-weighted mean plate degree.
    pub mean_degree: f64,
    /// Maximum degree among plates.
    pub max_degree: u32,
    /// Approximate triple-junction count (cells whose 1-ring spans ≥3 plates).
    pub triple_junctions: u32,
    /// Count of degree-1 plates (dangling, potential gaps/ends).
    pub degree_one: u32,
    /// Count of high-degree plates (≥5), potential overcrowded junctions.
    pub degree_ge5: u32,
}

/// Build plate adjacency from classified boundary edges.
/// Returns an undirected adjacency map from plate_id → neighbor set.
pub fn build_plate_adjacency(
    boundaries: &Boundaries,
    plates: &Plates,
) -> HashMap<u16, HashSet<u16>> {
    let mut adj: HashMap<u16, HashSet<u16>> = HashMap::new();
    for &(u, v, _cls) in &boundaries.edges {
        let pu = plates.plate_id[u as usize];
        let pv = plates.plate_id[v as usize];
        if pu == pv {
            continue;
        }
        adj.entry(pu).or_default().insert(pv);
        adj.entry(pv).or_default().insert(pu);
    }
    adj
}

/// Count cells whose 1-ring neighborhood spans at least 3 distinct plate labels.
pub fn count_triple_junction_like(grid: &Grid, plates: &Plates) -> u32 {
    let mut count = 0u32;
    for i in 0..grid.cells {
        let mut uniq: HashSet<u16> = HashSet::with_capacity(8);
        uniq.insert(plates.plate_id[i]);
        for &nj in &grid.n1[i] {
            uniq.insert(plates.plate_id[nj as usize]);
            if uniq.len() >= 3 {
                count += 1;
                break;
            }
        }
    }
    count
}

/// Compute high-level connectivity and junction stats for the current boundary/plate state.
pub fn compute_stats(grid: &Grid, plates: &Plates, boundaries: &Boundaries) -> PlateNetworkStats {
    // Plate set present
    let mut present: HashSet<u16> = HashSet::new();
    for &pid in &plates.plate_id {
        present.insert(pid);
    }
    let n_plates = present.len() as u32;

    // Adjacency and degrees
    let adj = build_plate_adjacency(boundaries, plates);
    let mut degree_sum = 0u64;
    let mut degree_max = 0u32;
    let mut edges_undirected = 0u64;
    let mut degree_one = 0u32;
    let mut degree_ge5 = 0u32;
    for (&pid, nbrs) in &adj {
        let deg = nbrs.len() as u32;
        degree_sum += deg as u64;
        if deg > degree_max {
            degree_max = deg;
        }
        // Count unique edges by summing degrees and dividing by 2 later
        edges_undirected += deg as u64;
        // Ensure isolated plates appear in adj map
        let _ = present.contains(&pid);
        if deg == 1 {
            degree_one += 1;
        }
        if deg >= 5 {
            degree_ge5 += 1;
        }
    }
    let edges = (edges_undirected / 2) as u32;
    let mean_degree = if n_plates > 0 { (degree_sum as f64) / (n_plates as f64) } else { 0.0 };

    // Connected components over plate graph
    let mut comps = 0u32;
    let mut visited: HashSet<u16> = HashSet::new();
    for &pid in &present {
        if visited.contains(&pid) {
            continue;
        }
        comps += 1;
        let mut q: VecDeque<u16> = VecDeque::new();
        visited.insert(pid);
        q.push_back(pid);
        while let Some(p) = q.pop_front() {
            if let Some(nbrs) = adj.get(&p) {
                for &qpid in nbrs {
                    if visited.insert(qpid) {
                        q.push_back(qpid);
                    }
                }
            }
        }
    }

    let tj = count_triple_junction_like(grid, plates);
    PlateNetworkStats {
        num_plates: n_plates,
        edges,
        components: comps,
        mean_degree,
        max_degree: degree_max,
        triple_junctions: tj,
        degree_one,
        degree_ge5,
    }
}

/// Compact one-line summary string for logs
pub fn summarize(stats: &PlateNetworkStats) -> String {
    format!(
        "plates={} edges={} comps={} deg_mean={:.2} deg_max={} deg1={} deg>=5={} tj~{}",
        stats.num_plates,
        stats.edges,
        stats.components,
        stats.mean_degree,
        stats.max_degree,
        stats.degree_one,
        stats.degree_ge5,
        stats.triple_junctions
    )
}

/// Vertex-degree health over the boundary graph (cells as vertices, boundary edges as edges).
#[derive(Debug, Clone, Copy, Default)]
pub struct BoundaryHealthStats {
    /// Number of boundary vertices with degree 1 (chain ends). Open chains ≈ deg1/2.
    pub deg1: u32,
    /// Number of boundary vertices with degree 2 (well-formed chains/loops interior).
    pub deg2: u32,
    /// Number of boundary vertices with degree ≥3 (junctions/crossings).
    pub deg3p: u32,
    /// Approximate number of open chains (= deg1/2, integer division).
    pub open_chains: u32,
}

/// Compute boundary vertex-degree histogram for health diagnostics.
pub fn compute_boundary_health(boundaries: &Boundaries) -> BoundaryHealthStats {
    let mut deg: HashMap<u32, u32> = HashMap::new();
    for &(u, v, _cls) in &boundaries.edges {
        *deg.entry(u).or_insert(0) += 1;
        *deg.entry(v).or_insert(0) += 1;
    }
    let mut d1 = 0u32;
    let mut d2 = 0u32;
    let mut d3p = 0u32;
    for (_cell, d) in deg {
        match d {
            0 => {}
            1 => d1 += 1,
            2 => d2 += 1,
            _ => d3p += 1,
        }
    }
    BoundaryHealthStats { deg1: d1, deg2: d2, deg3p: d3p, open_chains: d1 / 2 }
}

/// Compact one-line summary for boundary health.
pub fn summarize_boundary(h: &BoundaryHealthStats) -> String {
    format!(
        "boundary: deg1={} deg2={} deg3+={} open_chains~{}",
        h.deg1, h.deg2, h.deg3p, h.open_chains
    )
}
