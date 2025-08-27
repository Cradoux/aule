// Keep imports minimal to ease no-std migration if needed

pub const EPS_ROLLOVER: f32 = 1.0e-7;
pub const EPS_UPPER: f32 = 1.0e-6;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };
    #[must_use]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
    #[must_use]
    pub fn dot(self, o: Self) -> f32 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    #[must_use]
    pub fn length(self) -> f32 {
        self.dot(self).sqrt()
    }
    #[must_use]
    pub fn normalized(self) -> Self {
        let l = self.length();
        if l == 0.0 {
            self
        } else {
            Self::new(self.x / l, self.y / l, self.z / l)
        }
    }
}

use core::ops::{Add, Mul, Sub};

impl Add for Vec3 {
    type Output = Self;
    #[must_use]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[must_use]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f32> for Vec3 {
    type Output = Self;
    #[must_use]
    fn mul(self, rhs: f32) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

pub type FaceId = u32;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct FaceGeom {
    /// CCW from outside. Unit-length vertices on the sphere.
    pub a: Vec3,
    pub b: Vec3,
    pub c: Vec3,
    /// Outward unit normal of the plane (matches (b-a) × (c-a) after winding fix).
    pub n: Vec3,
    /// Neighbor ids opposite A, B, C (i.e., across BC, CA, AB).
    pub neighbor_opp: [FaceId; 3],
}

impl FaceGeom {
    #[inline]
    #[must_use]
    #[allow(clippy::similar_names)]
    pub fn ensure_outward_ccw(mut self) -> Self {
        // Outward iff dot(normal, centroid) > 0
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        let mut n = ab.cross(ac).normalized();
        let centroid = self.a + self.b + self.c;
        if n.dot(centroid) < 0.0 {
            // Flip winding B <-> C
            std::mem::swap(&mut self.b, &mut self.c);
            let edge_ab = self.b - self.a;
            let edge_ac = self.c - self.a;
            n = edge_ab.cross(edge_ac).normalized();
        }
        self.n = n;
        self
    }

    /// Compute gnomonic barycentrics for point p on the unit sphere using tangent basis at face centroid.
    /// Returns (α,β,γ). This reduces distortion near edges.
    #[must_use]
    pub fn barycentrics_gnomonic(&self, point: Vec3) -> [f32; 3] {
        // Face centroid (approx) and orthonormal basis in tangent plane
        let centroid_dir = (self.a + self.b + self.c).normalized();
        let tangent_u = (self.a - centroid_dir * self.a.dot(centroid_dir)).normalized();
        let tangent_v = centroid_dir.cross(tangent_u).normalized();
        // Project triangle corners to tangent plane using gnomonic projection
        let proj = |vertex: Vec3| -> (f32, f32) {
            let denom = vertex.dot(centroid_dir).max(1e-9);
            let inv_denom = 1.0 / denom;
            let proj_x = vertex.dot(tangent_u) * inv_denom;
            let proj_y = vertex.dot(tangent_v) * inv_denom;
            (proj_x, proj_y)
        };
        let (ax, ay) = proj(self.a);
        let (bx, by) = proj(self.b);
        let (cx, cy) = proj(self.c);
        let (px, py) = proj(point);
        // Planar barycentrics in tangent plane
        let edge0_x = bx - ax;
        let edge0_y = by - ay;
        let edge1_x = cx - ax;
        let edge1_y = cy - ay;
        let delta_x = px - ax;
        let delta_y = py - ay;
        let m_e0e0 = edge0_x * edge0_x + edge0_y * edge0_y;
        let m_e0e1 = edge0_x * edge1_x + edge0_y * edge1_y;
        let m_e1e1 = edge1_x * edge1_x + edge1_y * edge1_y;
        let m_pe0 = delta_x * edge0_x + delta_y * edge0_y;
        let m_pe1 = delta_x * edge1_x + delta_y * edge1_y;
        let denom = (m_e0e0 * m_e1e1 - m_e0e1 * m_e0e1).max(1e-12);
        let beta = (m_e1e1 * m_pe0 - m_e0e1 * m_pe1) / denom;
        let gamma = (m_e0e0 * m_pe1 - m_e0e1 * m_pe0) / denom;
        let alpha = 1.0 - beta - gamma;
        [alpha, beta, gamma]
    }
}

/// Argmax `dot(N_f, p)`. Assumes all `FaceGeom` have outward n.
/// # Panics
/// Panics if the best face index does not fit into `u32` (should never happen).
#[inline]
#[must_use]
pub fn pick_face(p: Vec3, faces: &[FaceGeom]) -> FaceId {
    let mut best_i: usize = 0;
    let mut best_d = f32::NEG_INFINITY;
    for (i, f) in faces.iter().enumerate() {
        let d = f.n.dot(p);
        if d > best_d || ((d - best_d).abs() <= f32::EPSILON && i < best_i) {
            best_d = d;
            best_i = i;
        }
    }
    core::convert::TryFrom::try_from(best_i).expect("face index fits in u32")
}

/// Project p onto face plane and compute planar barycentrics (α,β,γ) wrt (A,B,C).
#[inline]
#[must_use]
pub fn barycentrics_plane(p: Vec3, f: &FaceGeom) -> [f32; 3] {
    // Plane projection
    let to_plane = f.n * f.n.dot(p); // component along n
    let p_proj = p - to_plane;

    // Build local basis at A
    let v0 = f.b - f.a;
    let v1 = f.c - f.a;
    let v2 = p_proj - f.a;

    let dot00 = v0.dot(v0);
    let dot01 = v0.dot(v1);
    let dot02 = v0.dot(v2);
    let dot11 = v1.dot(v1);
    let dot12 = v1.dot(v2);

    let denom = dot00 * dot11 - dot01 * dot01;
    // If denom is tiny, fallback (degenerate triangle shouldn't happen, but guard anyway)
    if denom.abs() < 1e-20 {
        return [1.0, 0.0, 0.0];
    }

    // Solve v2 = u*v0 + v*v1
    let u = (dot11 * dot02 - dot01 * dot12) / denom;
    let v = (dot00 * dot12 - dot01 * dot02) / denom;

    // Barycentrics: α at A, β at B, γ at C
    let alpha = 1.0 - u - v;
    let beta = u;
    let gamma = v;
    [alpha, beta, gamma]
}

/// Return (face, w=[α,β,γ], kneg) after one-hop rollover if needed.
/// `eps_rollover` allows tiny negatives due to FP.
#[inline]
#[must_use]
pub fn one_hop_rollover(
    faces: &[FaceGeom],
    face: FaceId,
    w: [f32; 3],
    eps_rollover: f32,
) -> (FaceId, [f32; 3], u32) {
    let f_id = face;
    // Find most-negative component
    let (kneg, min_val) = {
        let mut idx = 0u32;
        let mut val = w[0];
        if w[1] < val {
            idx = 1;
            val = w[1];
        }
        if w[2] < val {
            idx = 2;
            val = w[2];
        }
        (idx, val)
    };

    if min_val < -eps_rollover {
        let next = faces[f_id as usize].neighbor_opp[kneg as usize];
        let f2 = &faces[next as usize];
        let w2 = barycentrics_plane(
            // Recompute on neighbor plane (use same 3D p recovered from barycentrics on f)
            // We cannot recover exact p here without it; assume caller recomputes from lon/lat.
            // Practical path: caller has p and passes it here; adjust signature if desired.
            // For now, we assume caller recomputes barycentrics on f2 immediately after.
            // This function just returns the intent (kneg & neighbor).
            // To keep it stateless, we return (next, [NaN;3]) to signal "recompute".
            // ---- Implementation note ----
            // In practice, we don’t want NaNs. So we keep the signature but expect caller to
            // immediately recompute barycentrics_plane(p, faces[next]).
            Vec3::ZERO, // placeholder, not used
            f2,
        );
        // We cannot compute w2 without p. Return marker:
        return (next, w2, kneg);
    }
    (f_id, w, kneg)
}

/// Lattice mapping: u ← α·F, v ← β·F; reflect if iu+iv>F−1; “upper” if
/// frac(u)+frac(v) ≥ `1−eps_upper`.
#[inline]
#[must_use]
pub fn lattice_from_ab(alpha: f32, beta: f32, f: u32, eps_upper: f32) -> (u32, u32, bool) {
    #[allow(clippy::cast_precision_loss)]
    let fu = alpha * (f as f32);
    #[allow(clippy::cast_precision_loss)]
    let fv = beta * (f as f32);

    let iu = fu.floor();
    let iv = fv.floor();

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let mut u = iu as u32;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let mut v = iv as u32;

    let mut ru = fu - iu;
    let mut rv = fv - iv;

    if u + v > f.saturating_sub(1) {
        // reflect across diagonal into principal triangle
        u = f.saturating_sub(1) - u;
        v = f.saturating_sub(1) - v;
        ru = 1.0 - ru;
        rv = 1.0 - rv;
    }

    // Strict '>' for upper; diagonal belongs to lower
    let upper = (ru + rv) > (1.0 - eps_upper);
    (u, v, upper)
}

/// Triangle indexing: MUST match the engine’s existing scheme exactly.
/// Replace the body with your engine’s current implementation.
pub struct TriIndex;

impl TriIndex {
    /// Compute triangle id within a face from (u,v,upper).
    /// Total triangles per face = 2*F*F.
    #[allow(unused_variables)]
    #[must_use]
    pub fn tri_index(f: u32, u: u32, v: u32, upper: bool) -> u32 {
        // Canonical per-face formula (matches WGSL):
        // tri_local = v*(2F - v) + 2*u + (upper?1:0)
        v * (2 * f - v) + 2 * u + u32::from(upper)
    }
}

/// Convert spherical (lon, lat) in radians to a unit vector on the sphere.
#[inline]
#[must_use]
pub fn sph_to_unit(lon: f32, lat: f32) -> Vec3 {
    let (slon, clon) = lon.sin_cos();
    let (slat, clat) = lat.sin_cos();
    // Canonical mapping (note the minus on z):
    // x = cos(lat) * cos(lon), y = sin(lat), z = -cos(lat) * sin(lon)
    Vec3 { x: clon * clat, y: slat, z: -slon * clat }
}
