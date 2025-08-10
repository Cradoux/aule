//! Geodesic icosphere grid and adjacency.

pub mod cache;
pub mod tile;

use smallvec::SmallVec;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::f64::consts::PI;

/// Grid data for an icosphere dual (Voronoi) discretization on the unit sphere.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    /// Number of cells (equal to number of original icosphere vertices)
    pub cells: usize,
    /// Cell center positions in XYZ (unit sphere)
    pub pos_xyz: Vec<[f32; 3]>,
    /// Cell center latitude/longitude in radians
    pub latlon: Vec<[f32; 2]>,
    /// Dual cell spherical area on unit sphere
    pub area: Vec<f32>,
    /// 1-ring neighbor cell indices (pentagon:5, hexagon:6)
    pub n1: Vec<SmallVec<[u32; 6]>>,
    /// 2-ring neighbor indices (unique, excludes self)
    pub n2: Vec<SmallVec<[u32; 12]>>,
    /// Icosphere subdivision frequency
    pub frequency: u32,
}

impl Grid {
    /// Build an icosphere dual grid with geodesic frequency F (Class-I, u=F, v=0).
    /// Geometry and area are computed in f64 then stored as f32.
    pub fn new(frequency: u32) -> Self {
        let f = frequency.max(1); // treat F=0 as F=1 (base icosahedron)

        // Canonical icosahedron (unit sphere)
        let (v0, faces) = icosahedron();

        // Canonical ID space built in two phases: assign IDs deterministically, then emit positions.
        // Vertex definition per ID for later position emission
        enum VertexDef {
            Corner(usize),
            Edge { lo: usize, hi: usize, t_from_lo: u32 },
            Interior { face: usize, i: u32, j: u32 },
        }
        let mut id_defs: Vec<VertexDef> = Vec::new();

        // Corner map: base corner index -> global id
        let mut corner_map: HashMap<usize, u32> = HashMap::new();
        for i in 0..v0.len() {
            let id = id_defs.len() as u32;
            corner_map.insert(i, id);
            id_defs.push(VertexDef::Corner(i));
        }

        // Edge map: (lo, hi, t_from_lo) with 1..F-1
        #[derive(Hash, Eq, PartialEq, Copy, Clone)]
        struct EdgeKey(u32, u32, u32);
        let mut edge_map: HashMap<EdgeKey, u32> = HashMap::new();

        // Utility to produce an EdgeKey from an ordered pair (a,b) and a parameter t from A
        let edge_key_from = |a_id: usize, b_id: usize, t_from_a: u32| -> (EdgeKey, u32) {
            let (lo_usize, hi_usize, t_from_lo) = if a_id < b_id {
                (a_id, b_id, t_from_a)
            } else {
                // distance from lo when lo is b
                (b_id, a_id, f - t_from_a)
            };
            (EdgeKey(lo_usize as u32, hi_usize as u32, t_from_lo), t_from_lo)
        };

        // ID assignment grid per face
        let mut tris: Vec<[u32; 3]> = Vec::new();
        for (fid, face) in faces.iter().enumerate() {
            let (a_id, b_id, c_id) = (face[0], face[1], face[2]);
            // Build grid_ids for this face
            let mut grid_ids: Vec<Vec<u32>> = Vec::with_capacity((f as usize) + 1);
            for i in 0..=f {
                let mut row: Vec<u32> = Vec::with_capacity((f - i) as usize + 1);
                for j in 0..=(f - i) {
                    let k = f - i - j;
                    let id = if j == 0 && k == 0 {
                        corner_map[&a_id]
                    } else if i == 0 && k == 0 {
                        corner_map[&b_id]
                    } else if i == 0 && j == 0 {
                        corner_map[&c_id]
                    } else if k == 0 {
                        // AB edge, t from A
                        let t = j;
                        let (ekey, t_from_lo) = edge_key_from(a_id, b_id, t);
                        if let Some(&eid) = edge_map.get(&ekey) {
                            eid
                        } else {
                            let eid = id_defs.len() as u32;
                            edge_map.insert(ekey, eid);
                            id_defs.push(VertexDef::Edge {
                                lo: ekey.0 as usize,
                                hi: ekey.1 as usize,
                                t_from_lo,
                            });
                            eid
                        }
                    } else if j == 0 {
                        // C-A edge, local t = i (distance from C toward A)
                        let t = i;
                        let (ekey, t_from_lo) = edge_key_from(c_id, a_id, t);
                        if let Some(&eid) = edge_map.get(&ekey) {
                            eid
                        } else {
                            let eid = id_defs.len() as u32;
                            edge_map.insert(ekey, eid);
                            id_defs.push(VertexDef::Edge {
                                lo: ekey.0 as usize,
                                hi: ekey.1 as usize,
                                t_from_lo,
                            });
                            eid
                        }
                    } else if i == 0 {
                        // B-C edge, t from B
                        let t = k;
                        let (ekey, t_from_lo) = edge_key_from(b_id, c_id, t);
                        if let Some(&eid) = edge_map.get(&ekey) {
                            eid
                        } else {
                            let eid = id_defs.len() as u32;
                            edge_map.insert(ekey, eid);
                            id_defs.push(VertexDef::Edge {
                                lo: ekey.0 as usize,
                                hi: ekey.1 as usize,
                                t_from_lo,
                            });
                            eid
                        }
                    } else {
                        // Interior unique to this face
                        let nid = id_defs.len() as u32;
                        id_defs.push(VertexDef::Interior { face: fid, i, j });
                        nid
                    };
                    row.push(id);
                }
                grid_ids.push(row);
            }

            // Canonical two-loop triangulation: exactly F^2 triangles per face (no guards)
            let f_usize = f as usize;
            let idx = |ii: usize, jj: usize| -> u32 { grid_ids[ii][jj] };
            #[cfg(debug_assertions)]
            let mut face_tri_count: usize = 0;
            #[cfg(debug_assertions)]
            let mut face_area_sum: f64 = 0.0;
            #[cfg(debug_assertions)]
            let get_pos = |vid: usize| -> [f64; 3] {
                match id_defs[vid] {
                    VertexDef::Corner(c) => v0[c],
                    VertexDef::Edge { lo, hi, t_from_lo } => {
                        let lo_v = v0[lo];
                        let hi_v = v0[hi];
                        let tt = t_from_lo as f64 / f as f64;
                        normalize3([
                            lo_v[0] * (1.0 - tt) + hi_v[0] * tt,
                            lo_v[1] * (1.0 - tt) + hi_v[1] * tt,
                            lo_v[2] * (1.0 - tt) + hi_v[2] * tt,
                        ])
                    }
                    VertexDef::Interior { face, i, j } => {
                        let (aa, bb, cc) =
                            (v0[faces[face][0]], v0[faces[face][1]], v0[faces[face][2]]);
                        let k = f - i - j;
                        let fi = f as f64;
                        let w_a = i as f64 / fi;
                        let w_b = j as f64 / fi;
                        let w_c = k as f64 / fi;
                        normalize3([
                            aa[0] * w_a + bb[0] * w_b + cc[0] * w_c,
                            aa[1] * w_a + bb[1] * w_b + cc[1] * w_c,
                            aa[2] * w_a + bb[2] * w_b + cc[2] * w_c,
                        ])
                    }
                }
            };
            // DOWN triangles: [i,j],[i+1,j],[i,j+1]
            for i in 0..f_usize {
                for j in 0..(f_usize - i) {
                    let a = idx(i, j);
                    let b = idx(i + 1, j);
                    let c = idx(i, j + 1);
                    #[cfg(debug_assertions)]
                    {
                        let pa = get_pos(a as usize);
                        let pb = get_pos(b as usize);
                        let pc = get_pos(c as usize);
                        face_area_sum += spherical_triangle_area(pa, pb, pc);
                        face_tri_count += 1;
                    }
                    tris.push([a, b, c]);
                }
            }
            // UP triangles: [i+1,j],[i+1,j+1],[i,j+1]
            for i in 0..f_usize {
                for j in 0..(f_usize - i - 1) {
                    let b = idx(i + 1, j);
                    let d = idx(i + 1, j + 1);
                    let c = idx(i, j + 1);
                    #[cfg(debug_assertions)]
                    {
                        let pb = get_pos(b as usize);
                        let pd = get_pos(d as usize);
                        let pc = get_pos(c as usize);
                        face_area_sum += spherical_triangle_area(pb, pd, pc);
                        face_tri_count += 1;
                    }
                    tris.push([b, d, c]);
                }
            }
            #[cfg(debug_assertions)]
            {
                debug_assert_eq!(face_tri_count as u32, f * f);
                let a = v0[a_id];
                let b = v0[b_id];
                let c = v0[c_id];
                let face_true = spherical_triangle_area(a, b, c);
                let rel = (face_area_sum - face_true).abs() / face_true.max(1e-18);
                debug_assert!(rel < 1e-10, "face {} area rel err {}", fid, rel);
            }
        }

        // Emit positions for all IDs
        let n_vertices = id_defs.len();
        let mut global_pos: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0]; n_vertices];
        for (id, def) in id_defs.iter().enumerate() {
            let pos = match *def {
                VertexDef::Corner(c) => v0[c],
                VertexDef::Edge { lo, hi, t_from_lo } => {
                    let lo_v = v0[lo];
                    let hi_v = v0[hi];
                    let tt = t_from_lo as f64 / f as f64;
                    normalize3([
                        lo_v[0] * (1.0 - tt) + hi_v[0] * tt,
                        lo_v[1] * (1.0 - tt) + hi_v[1] * tt,
                        lo_v[2] * (1.0 - tt) + hi_v[2] * tt,
                    ])
                }
                VertexDef::Interior { face, i, j } => {
                    let (a_id, b_id, c_id) = (faces[face][0], faces[face][1], faces[face][2]);
                    let (a, b, c) = (v0[a_id], v0[b_id], v0[c_id]);
                    let k = f - i - j;
                    let fi = f as f64;
                    let w_a = i as f64 / fi;
                    let w_b = j as f64 / fi;
                    let w_c = k as f64 / fi;
                    normalize3([
                        a[0] * w_a + b[0] * w_b + c[0] * w_c,
                        a[1] * w_a + b[1] * w_b + c[1] * w_c,
                        a[2] * w_a + b[2] * w_b + c[2] * w_c,
                    ])
                }
            };
            global_pos[id] = pos;
        }

        // Build 1-ring adjacency from triangles
        let n_vertices = global_pos.len();
        let mut n1_sets: Vec<HashSet<u32>> = vec![HashSet::new(); n_vertices];
        for t in &tris {
            let a = t[0] as usize;
            let b = t[1] as usize;
            let c = t[2] as usize;
            n1_sets[a].insert(t[1]);
            n1_sets[a].insert(t[2]);
            n1_sets[b].insert(t[0]);
            n1_sets[b].insert(t[2]);
            n1_sets[c].insert(t[0]);
            n1_sets[c].insert(t[1]);
        }

        // Build 2-ring as unique neighbors of neighbors excluding self
        let mut n1: Vec<SmallVec<[u32; 6]>> = Vec::with_capacity(n_vertices);
        let mut n2: Vec<SmallVec<[u32; 12]>> = Vec::with_capacity(n_vertices);
        for vid in 0..n_vertices {
            let mut n1_sorted: Vec<u32> = n1_sets[vid].iter().copied().collect();
            n1_sorted.sort_unstable();
            let mut n2_set: BTreeSet<u32> = BTreeSet::new();
            for &u in &n1_sorted {
                for &v in &n1_sets[u as usize] {
                    if v != vid as u32 && !n1_sets[vid].contains(&v) {
                        n2_set.insert(v);
                    }
                }
            }
            let mut n1_sv: SmallVec<[u32; 6]> = SmallVec::new();
            for x in n1_sorted {
                n1_sv.push(x);
            }
            let mut n2_sv: SmallVec<[u32; 12]> = SmallVec::new();
            for x in n2_set {
                n2_sv.push(x);
            }
            n1.push(n1_sv);
            n2.push(n2_sv);
        }

        // Areas: distribute each spherical triangle's area equally to its three vertices
        let mut area_acc: Vec<f64> = vec![0.0; n_vertices];
        for t in &tris {
            let a = global_pos[t[0] as usize];
            let b = global_pos[t[1] as usize];
            let c = global_pos[t[2] as usize];
            let tri_area = spherical_triangle_area(a, b, c);
            let share = tri_area / 3.0;
            area_acc[t[0] as usize] += share;
            area_acc[t[1] as usize] += share;
            area_acc[t[2] as usize] += share;
        }

        // Validate partitioning of sphere area ~ 4π
        let sum_area: f64 = area_acc.iter().sum();
        let rel_err = (sum_area - 4.0 * PI).abs() / (4.0 * PI);
        debug_assert!(rel_err < 1e-6, "area partition error: {rel_err}");

        // Convert to output types
        let cells = n_vertices;
        let mut pos_xyz: Vec<[f32; 3]> = Vec::with_capacity(cells);
        let mut latlon: Vec<[f32; 2]> = Vec::with_capacity(cells);
        let mut area: Vec<f32> = Vec::with_capacity(cells);
        for p in &global_pos {
            pos_xyz.push([p[0] as f32, p[1] as f32, p[2] as f32]);
            let ll = xyz_to_latlon_f64(*p);
            latlon.push([ll[0] as f32, ll[1] as f32]);
        }
        for a in &area_acc {
            area.push(*a as f32);
        }

        Self { cells, pos_xyz, latlon, area, n1, n2, frequency }
    }
}

/// Utility: convert Cartesian to lat/lon (radians) using f64
fn xyz_to_latlon_f64(p: [f64; 3]) -> [f64; 2] {
    let (x, y, z) = (p[0], p[1], p[2]);
    let lat = z.atan2((x * x + y * y).sqrt());
    let lon = y.atan2(x);
    [lat, lon]
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn triple(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    dot(cross(a, b), c)
}

/// Spherical triangle area on unit sphere using the robust vector formula.
fn spherical_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let numerator = triple(a, b, c).abs();
    let denom = 1.0 + dot(a, b) + dot(b, c) + dot(c, a);
    2.0 * (numerator.atan2(denom))
}

/// Canonical icosahedron vertices (unit sphere) and faces (CCW)
fn icosahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let phi = (1.0 + 5.0_f64.sqrt()) * 0.5;
    let a = 1.0;
    let b = 1.0 / phi;
    let verts = vec![
        normalize3([-a, b, 0.0]),
        normalize3([a, b, 0.0]),
        normalize3([-a, -b, 0.0]),
        normalize3([a, -b, 0.0]),
        normalize3([0.0, -a, b]),
        normalize3([0.0, a, b]),
        normalize3([0.0, -a, -b]),
        normalize3([0.0, a, -b]),
        normalize3([b, 0.0, -a]),
        normalize3([b, 0.0, a]),
        normalize3([-b, 0.0, -a]),
        normalize3([-b, 0.0, a]),
    ];
    // Faces (indices into verts). Canonical CCW set.
    let faces = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    (verts, faces)
}
