#![forbid(unsafe_code)]
#![deny(clippy::all, clippy::pedantic)]

pub mod icosa;
mod math;

pub use icosa::{build_face_table, GpuFaceGeom, GpuNeighbors};
pub use math::{
    barycentrics_plane, lattice_from_ab, one_hop_rollover, pick_face, sph_to_unit, FaceGeom,
    FaceId, TriIndex, Vec3, EPS_ROLLOVER, EPS_UPPER,
};

/// Choose f32 everywhere to mirror WGSL exactly.
pub const FLOAT: &str = "f32";

/// Quantize α,β to lattice (u,v) and compute tri id using the engine’s row indexing.
/// IMPORTANT: replace the body of `TriIndex::tri_index` in math.rs with your engine’s formula.
/// This function is just a convenience façade that calls the shared implementation.
#[must_use]
pub fn quantize_and_index(alpha: f32, beta: f32, f: u32, eps_upper: f32) -> (u32, u32, bool, u32) {
    let (u, v, upper) = lattice_from_ab(alpha, beta, f, eps_upper);
    let tri = TriIndex::tri_index(f, u, v, upper);
    (u, v, upper, tri)
}
