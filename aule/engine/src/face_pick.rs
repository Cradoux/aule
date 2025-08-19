//! CPU face/triangle selection using the shared aule-geo crate.
//! Safe to call from parity/debug paths without affecting lattice construction.

use aule_geo::{
    barycentrics_plane, build_face_table, lattice_from_ab, one_hop_rollover, pick_face, sph_to_unit as geo_sph_to_unit,
    FaceGeom, FaceId, TriIndex, Vec3 as GeoVec3, EPS_ROLLOVER, EPS_UPPER,
};

/// Minimal local Vec3 shim so the engine doesn’t pull in any extra math deps.
#[derive(Clone, Copy, Debug)]
pub struct Vec3 {
    /// X component
    pub x: f32,
    /// Y component
    pub y: f32,
    /// Z component
    pub z: f32,
}
impl Vec3 {
    /// Construct from components.
    #[inline]
    pub fn new(x: f32, y: f32, z: f32) -> Self { Self { x, y, z } }
    /// Dot product with another vector.
    #[inline]
    pub fn dot(self, o: Self) -> f32 { self.x * o.x + self.y * o.y + self.z * o.z }
    /// Euclidean length.
    #[inline]
    pub fn len(self) -> f32 { self.dot(self).sqrt() }
    /// Unit-normalize, returning the input if zero length.
    #[inline]
    pub fn norm(self) -> Self {
        let l = self.len();
        if l == 0.0 { self } else { Self::new(self.x / l, self.y / l, self.z / l) }
    }
}

/// Convert lon/lat (radians) to unit vector (ECEF). Uses shared implementation.
#[inline]
pub fn sph_to_unit(lon: f32, lat: f32) -> Vec3 {
    let v = geo_sph_to_unit(lon, lat);
    Vec3::new(v.x, v.y, v.z)
}

/// Result of CPU pick.
#[derive(Clone, Copy, Debug)]
pub struct FaceTri {
    /// Face id 0..19
    pub face: FaceId,
    /// Triangle index within face (0..2*F*F-1)
    pub tri: u32,
    /// Lattice u coordinate
    pub u: u32,
    /// Lattice v coordinate
    pub v: u32,
    /// True if upper triangle
    pub upper: bool,
    /// Barycentrics on the final face
    pub w: [f32; 3],
    /// Which component was most-negative before rollover
    pub kneg: u32,
}

/// Stateless picker that owns the shared face table.
pub struct GeoPicker {
    faces: Vec<FaceGeom>,
}
impl GeoPicker {
    /// Create with canonical faces from `aule-geo`.
    pub fn new() -> Self { Self { faces: build_face_table() } }

    /// Access the face table.
    #[inline]
    pub fn faces(&self) -> &[FaceGeom] { &self.faces }

    /// Pick (face,tri) for a unit vector `p` at lattice resolution F.
    pub fn pick_face_tri_from_p(&self, p: Vec3, f_subdiv: u32) -> FaceTri {
        // 1) initial face by argmax dot(n, p)
        let p_geo = GeoVec3 { x: p.x, y: p.y, z: p.z };
        let f0 = pick_face(p_geo, &self.faces);
        // 2) planar barycentrics on that face
        let w0 = barycentrics_plane(p_geo, &self.faces[f0 as usize]);
        // 3) one-hop rollover if a component is < -EPS
        let (f1, w1, kneg) = {
            let (next, _w_marker, k) = one_hop_rollover(&self.faces, f0, w0, EPS_ROLLOVER);
            if next != f0 {
                let w2 = barycentrics_plane(p_geo, &self.faces[next as usize]);
                (next, w2, k)
            } else {
                (f0, w0, k)
            }
        };
        // 4) lattice quantization from α,β and tri indexing (matches WGSL)
        let (u, v, upper) = lattice_from_ab(w1[0], w1[1], f_subdiv, EPS_UPPER);
        let tri = TriIndex::tri_index(f_subdiv, u, v, upper);
        FaceTri { face: f1, tri, u, v, upper, w: w1, kneg }
    }

    /// Convenience: pick from lon/lat in radians.
    pub fn pick_face_tri_from_lonlat(&self, lon: f32, lat: f32, f_subdiv: u32) -> FaceTri {
        let p = sph_to_unit(lon, lat).norm();
        self.pick_face_tri_from_p(p, f_subdiv)
    }
}


