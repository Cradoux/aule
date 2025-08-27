//! 3D globe rendering: shared-vertex geodesic mesh and pipeline.

pub mod mesh;
pub mod orbit_cam;
pub mod pipeline;

pub use mesh::{build_globe_mesh, GlobeMesh};
pub use orbit_cam::OrbitCamera;
pub use pipeline::GlobeRenderer;
