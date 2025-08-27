//! Viewer-side CPU face/tri picker using the shared aule_geo crate.
//! Used by parity and "Force CPU face pick (debug)".

use aule_geo::{
    build_face_table, lattice_from_ab, pick_face, sph_to_unit as sph_to_unit_geo, FaceGeom, FaceId,
    TriIndex, Vec3 as GeoVec3, EPS_ROLLOVER, EPS_UPPER,
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
        // Use gnomonic barycentrics for robustness at edges
        let w0 = self.faces[f0 as usize].barycentrics_gnomonic(p3);

        // one-hop rollover only when exactly one component is within EPS of 0 and is the minimum
        let (f1, w1, kneg) = {
            let mut idx = 0usize;
            let mut vmin = w0[0];
            if w0[1] < vmin {
                vmin = w0[1];
                idx = 1;
            }
            if w0[2] < vmin {
                idx = 2;
            }
            if vmin.abs() <= EPS_ROLLOVER
                && (w0[(idx + 1) % 3] > EPS_ROLLOVER)
                && (w0[(idx + 2) % 3] > EPS_ROLLOVER)
            {
                let next = self.faces[f0 as usize].neighbor_opp[idx] as u32;
                let w2 = self.faces[next as usize].barycentrics_gnomonic(p3);
                (next, w2, idx as u32)
            } else {
                (f0, w0, idx as u32)
            }
        };

        let (u, v, upper) = lattice_from_ab(w1[0], w1[1], f_subdiv, EPS_UPPER);
        let tri = TriIndex::tri_index(f_subdiv, u, v, upper);
        FaceTri { face: f1, tri, u, v, upper, w: w1, kneg }
    }
}

impl Default for GeoPicker {
    fn default() -> Self {
        Self::new()
    }
}
