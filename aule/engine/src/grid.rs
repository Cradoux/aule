//! Geodesic icosphere grid and adjacency.

pub mod cache;

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

        // Vertex indexing maps for determinism
        let mut global_pos: Vec<[f64; 3]> = Vec::new();
        // Map base corner id -> global id
        let mut corner_map: HashMap<usize, u32> = HashMap::new();
        for i in 0..v0.len() {
            let id = global_pos.len() as u32;
            corner_map.insert(i, id);
            global_pos.push(v0[i]);
        }

        // Edge map: (min_corner, max_corner, t) with 1..F-1
        #[derive(Hash, Eq, PartialEq)]
        struct EdgeKey(u32, u32, u32);
        let mut edge_map: HashMap<EdgeKey, u32> = HashMap::new();

        // Face interior map: (face_id, i, j) with i>=1, j>=1, k>=1, i+j+k=F
        #[derive(Hash, Eq, PartialEq)]
        struct FaceKey(u32, u32, u32);
        let mut face_map: HashMap<FaceKey, u32> = HashMap::new();

        // Helper to add or fetch subdivided vertex on an edge
        let mut edge_point = |a_id: usize, b_id: usize, t: u32| -> u32 {
            let (lo, hi) =
                if a_id < b_id { (a_id as u32, b_id as u32) } else { (b_id as u32, a_id as u32) };
            let key = EdgeKey(lo, hi, t);
            if let Some(&gid) = edge_map.get(&key) {
                return gid;
            }
            let a = v0[a_id];
            let b = v0[b_id];
            let tt = t as f64 / f as f64;
            let p = normalize3([
                a[0] * (1.0 - tt) + b[0] * tt,
                a[1] * (1.0 - tt) + b[1] * tt,
                a[2] * (1.0 - tt) + b[2] * tt,
            ]);
            let gid = global_pos.len() as u32;
            edge_map.insert(key, gid);
            global_pos.push(p);
            gid
        };

        // Helper to add or fetch interior vertex
        let mut face_interior = |face_id: usize,
                                 i: u32,
                                 j: u32,
                                 k: u32,
                                 a: [f64; 3],
                                 b: [f64; 3],
                                 c: [f64; 3]|
         -> u32 {
            let key = FaceKey(face_id as u32, i, j);
            if let Some(&gid) = face_map.get(&key) {
                return gid;
            }
            let fi = f as f64;
            let w_a = i as f64 / fi;
            let w_b = j as f64 / fi;
            let w_c = k as f64 / fi;
            let p = normalize3([
                a[0] * w_a + b[0] * w_b + c[0] * w_c,
                a[1] * w_a + b[1] * w_b + c[1] * w_c,
                a[2] * w_a + b[2] * w_b + c[2] * w_c,
            ]);
            let gid = global_pos.len() as u32;
            face_map.insert(key, gid);
            global_pos.push(p);
            gid
        };

        // Build subdivided faces (small triangles) referencing global vertex ids
        let mut tris: Vec<[u32; 3]> = Vec::new();
        for (fid, face) in faces.iter().enumerate() {
            let (a_id, b_id, c_id) = (face[0], face[1], face[2]);
            let (a, b, c) = (v0[a_id], v0[b_id], v0[c_id]);

            // Generate grid points on this face
            // Barycentric i,j,k with i+j+k=F
            // Corners
            let a_gid = corner_map[&a_id];
            let b_gid = corner_map[&b_id];
            let c_gid = corner_map[&c_id];

            // Edge points
            // AB: t = 1..F-1
            let mut ab: Vec<u32> = Vec::new();
            for t in 1..f {
                ab.push(edge_point(a_id, b_id, t));
            }
            // BC
            let mut bc: Vec<u32> = Vec::new();
            for t in 1..f {
                bc.push(edge_point(b_id, c_id, t));
            }
            // CA
            let mut ca: Vec<u32> = Vec::new();
            for t in 1..f {
                ca.push(edge_point(c_id, a_id, t));
            }

            // Interior points per (i,j) with i>=1, j>=1, k>=1
            let mut interior_rows: Vec<Vec<u32>> = Vec::new();
            for i in 1..f {
                // i from 1..F-1
                let max_j = f - i;
                let mut row: Vec<u32> = Vec::new();
                for j in 1..max_j {
                    // j from 1..F-i-1, but include up to max_j-1 ensures k>=1
                    let k = f - i - j;
                    if k == 0 {
                        continue;
                    }
                    row.push(face_interior(fid, i, j, k, a, b, c));
                }
                if !row.is_empty() {
                    interior_rows.push(row);
                }
            }

            // Build small triangles in lexicographic order
            // We create a 2D grid layering from A->B along AB and A->C along AC.
            // Helper to index points on lattice rows
            // For each strip s = 0..F-1, we form two triangle types.

            // Precompute border arrays including corners for convenience
            let ab_full = {
                let mut v = Vec::with_capacity((f + 1) as usize);
                v.push(a_gid);
                v.extend(ab.iter().copied());
                v.push(b_gid);
                v
            };
            let bc_full = {
                let mut v = Vec::with_capacity((f + 1) as usize);
                v.push(b_gid);
                v.extend(bc.iter().copied());
                v.push(c_gid);
                v
            };
            let ca_full = {
                let mut v = Vec::with_capacity((f + 1) as usize);
                v.push(c_gid);
                v.extend(ca.iter().copied());
                v.push(a_gid);
                v
            };

            // Function to get point at barycentric (i,j,k) on this face mapped to global id
            let get_point = |i: u32, j: u32, k: u32| -> u32 {
                // corners
                if j == 0 && k == 0 {
                    return a_gid;
                }
                if i == 0 && k == 0 {
                    return b_gid;
                }
                if i == 0 && j == 0 {
                    return c_gid;
                }
                // edges
                if k == 0 {
                    // on AB: j in 1..F-1, i in 1..F-1 where j = t, i = F - t
                    let t = j; // distance from A
                    return ab_full[t as usize];
                }
                if j == 0 {
                    // on CA from C to A: k=t => index along CA_full from C
                    let t = k; // distance from C
                    return ca_full[t as usize];
                }
                if i == 0 {
                    // on BC from B to C: k=t => index along BC_full from B
                    let t = k; // distance from B to C
                    return bc_full[t as usize];
                }
                // interior: lookup from face_map using earlier insertion
                // We inserted only when i,j,k>=1
                let key = FaceKey(fid as u32, i, j);
                face_map[&key]
            };

            // Now create small triangles. For row r from 0..F-1, and column s from 0..(F-1-r)
            for i_b in 0..f {
                // i_b = barycentric i on A side for the strip
                let max_j = f - 1 - i_b;
                for j_b in 0..=max_j {
                    // Two triangles within the rhombus defined by (i,j) and (i+1,j)
                    if j_b < max_j {
                        // lower-left triangle: (i,j,k) - (i+1,j,k-1) - (i,j+1,k-1)
                        let i0 = i_b;
                        let j0 = j_b;
                        let k0 = f - i0 - j0;
                        if k0 == 0 {
                            continue;
                        }
                        let i1 = i_b + 1;
                        let j1 = j_b;
                        let k1 = f - i1 - j1;
                        let i2 = i_b;
                        let j2 = j_b + 1;
                        let k2 = f - i2 - j2;
                        if k1 == 0 || k2 == 0 {
                            continue;
                        }
                        let a0 = get_point(i0, j0, k0);
                        let b0 = get_point(i1, j1, k1);
                        let c0 = get_point(i2, j2, k2);
                        tris.push([a0, b0, c0]);
                    }
                    if i_b < f - 1 && j_b < max_j + 1 {
                        // upper-right triangle: (i+1,j,k-1) - (i+1,j+1,k-2) - (i,j+1,k-1)
                        if j_b < max_j {
                            let i1 = i_b + 1;
                            let j1 = j_b;
                            let k1 = f - i1 - j1;
                            let i2 = i_b + 1;
                            let j2 = j_b + 1;
                            let k2 = f - i2 - j2;
                            let i3 = i_b;
                            let j3 = j_b + 1;
                            let k3 = f - i3 - j3;
                            if k1 == 0 || k2 == 0 || k3 == 0 {
                                continue;
                            }
                            let a1 = get_point(i1, j1, k1);
                            let b1 = get_point(i2, j2, k2);
                            let c1 = get_point(i3, j3, k3);
                            tris.push([a1, b1, c1]);
                        }
                    }
                }
            }
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
        assert!(rel_err < 1e-6, "area partition error: {rel_err}");

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
    let mut verts = vec![
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
    (verts.drain(..).collect(), faces)
}
