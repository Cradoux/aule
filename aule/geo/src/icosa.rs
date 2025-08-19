use super::math::{FaceGeom, Vec3};
use std::collections::HashMap;

#[derive(Hash, Eq, PartialEq, Clone, Copy)]
struct Edge(usize, usize); // sorted

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GpuFaceGeom {
    pub a: [f32; 4],
    pub b: [f32; 4],
    pub c: [f32; 4],
    pub n: [f32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct GpuNeighbors {
    pub opp: [u32; 3], // [oppA, oppB, oppC]
}

/// Build canonical icosahedron faces with outward CCW and neighbor table.
///
/// # Panics
/// Panics if the neighbor detection fails to find exactly two faces sharing an edge.
#[must_use]
pub fn build_face_table() -> Vec<FaceGeom> {
    // Canonical icosahedron vertices (golden ratio φ)
    let phi = (1.0 + 5.0_f32.sqrt()) * 0.5;
    let inv_len = |v: Vec3| v.normalized();

    // 12 vertices (before normalization)
    let mut verts = vec![
        Vec3::new(-1.0, phi, 0.0),
        Vec3::new(1.0, phi, 0.0),
        Vec3::new(-1.0, -phi, 0.0),
        Vec3::new(1.0, -phi, 0.0),
        Vec3::new(0.0, -1.0, phi),
        Vec3::new(0.0, 1.0, phi),
        Vec3::new(0.0, -1.0, -phi),
        Vec3::new(0.0, 1.0, -phi),
        Vec3::new(phi, 0.0, -1.0),
        Vec3::new(phi, 0.0, 1.0),
        Vec3::new(-phi, 0.0, -1.0),
        Vec3::new(-phi, 0.0, 1.0),
    ];
    for v in &mut verts {
        *v = inv_len(*v);
    }

    // 20 faces (vertex indices), canonical ordering
    // Source: common icosahedron layout (matches many libs; we enforce outward CCW below).
    let faces_idx: [[usize; 3]; 20] = [
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

    let mut faces: Vec<FaceGeom> = faces_idx
        .iter()
        .map(|tri| {
            let f = FaceGeom {
                a: verts[tri[0]],
                b: verts[tri[1]],
                c: verts[tri[2]],
                n: Vec3::ZERO,
                neighbor_opp: [u32::MAX; 3],
            };
            f.ensure_outward_ccw()
        })
        .collect();

    // Build neighbor map via edge -> faces
    let mut edge_to_face: HashMap<Edge, Vec<usize>> = HashMap::new();

    for (fi, tri) in faces_idx.iter().enumerate() {
        let e0 = Edge(tri[1].min(tri[2]), tri[1].max(tri[2])); // edge BC (opp A)
        let e1 = Edge(tri[2].min(tri[0]), tri[2].max(tri[0])); // edge CA (opp B)
        let e2 = Edge(tri[0].min(tri[1]), tri[0].max(tri[1])); // edge AB (opp C)
        edge_to_face.entry(e0).or_default().push(fi);
        edge_to_face.entry(e1).or_default().push(fi);
        edge_to_face.entry(e2).or_default().push(fi);
    }

    // Fill neighbor_opp = [oppA(BC), oppB(CA), oppC(AB)]
    for (fi, tri) in faces_idx.iter().enumerate() {
        let look = |a: usize, b: usize| -> u32 {
            let e = Edge(a.min(b), a.max(b));
            let v = &edge_to_face[&e];
            let other = if v[0] == fi { v[1] } else { v[0] };
            u32::try_from(other).expect("face index fits in u32")
        };
        // Opposite A -> edge BC
        faces[fi].neighbor_opp[0] = look(tri[1], tri[2]);
        // Opposite B -> edge CA
        faces[fi].neighbor_opp[1] = look(tri[2], tri[0]);
        // Opposite C -> edge AB
        faces[fi].neighbor_opp[2] = look(tri[0], tri[1]);
    }

    faces
}

/// Convert CPU `FaceGeom` to GPU POD (std430-safe).
#[must_use]
pub fn to_gpu_faces(faces: &[FaceGeom]) -> (Vec<GpuFaceGeom>, Vec<GpuNeighbors>) {
    let pod_faces = faces
        .iter()
        .map(|f| GpuFaceGeom {
            a: [f.a.x, f.a.y, f.a.z, 0.0],
            b: [f.b.x, f.b.y, f.b.z, 0.0],
            c: [f.c.x, f.c.y, f.c.z, 0.0],
            n: [f.n.x, f.n.y, f.n.z, 0.0],
        })
        .collect::<Vec<_>>();
    let neigh = faces.iter().map(|f| GpuNeighbors { opp: f.neighbor_opp }).collect::<Vec<_>>();
    (pod_faces, neigh)
}
