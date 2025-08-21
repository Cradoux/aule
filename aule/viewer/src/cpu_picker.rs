//! Viewer-side CPU face/tri picker using the shared aule_geo crate.
//! Used by parity and "Force CPU face pick (debug)".

use aule_geo::{
    barycentrics_plane, build_face_table, lattice_from_ab, one_hop_rollover, pick_face,
    sph_to_unit as sph_to_unit_geo, FaceGeom, FaceId, TriIndex, Vec3 as GeoVec3, EPS_ROLLOVER,
    EPS_UPPER,
};

#[derive(Clone, Copy, Debug)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}
impl Vec3 {
    #[inline]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
    #[inline]
    pub fn norm(self) -> Self {
        let l = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if l == 0.0 {
            self
        } else {
            Self::new(self.x / l, self.y / l, self.z / l)
        }
    }
}

#[inline]
pub fn sph_to_unit(lon: f32, lat: f32) -> Vec3 {
    let v = sph_to_unit_geo(lon, lat);
    Vec3::new(v.x, v.y, v.z)
}

#[derive(Clone, Copy, Debug)]
pub struct FaceTri {
    pub face: FaceId,
    pub tri: u32,
    pub u: u32,
    pub v: u32,
    pub upper: bool,
    pub w: [f32; 3],
    pub kneg: u32,
}

pub struct GeoPicker {
    faces: Vec<FaceGeom>,
}
impl GeoPicker {
    pub fn new() -> Self {
        Self { faces: build_face_table() }
    }
    #[inline]
    pub fn faces(&self) -> &[FaceGeom] {
        &self.faces
    }

    pub fn pick_from_unit(&self, p: Vec3, f_subdiv: u32) -> FaceTri {
        let p3 = GeoVec3 { x: p.x, y: p.y, z: p.z };
        let f0 = pick_face(p3, &self.faces);
        let w0 = barycentrics_plane(p3, &self.faces[f0 as usize]);

        // one-hop rollover (shared signature without p → recompute if hopped)
        let (f1, w1, kneg) = {
            let (next, _marker, k) = one_hop_rollover(&self.faces, f0, w0, EPS_ROLLOVER);
            if next != f0 {
                let w2 = barycentrics_plane(p3, &self.faces[next as usize]);
                (next, w2, k)
            } else {
                (f0, w0, k)
            }
        };

        let (u, v, upper) = lattice_from_ab(w1[0], w1[1], f_subdiv, EPS_UPPER);
        let tri = TriIndex::tri_index(f_subdiv, u, v, upper);
        FaceTri { face: f1, tri, u, v, upper, w: w1, kneg }
    }
}
