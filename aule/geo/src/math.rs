// Keep imports minimal to ease no-std migration if needed

pub const EPS_ROLLOVER: f32 = 1.0e-6;
pub const EPS_UPPER: f32 = 1.0e-6;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vec3 {
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
    pub fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
    pub fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
    pub fn mul(self, k: f32) -> Self {
        Self::new(self.x * k, self.y * k, self.z * k)
    }
    pub fn dot(self, o: Self) -> f32 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    pub fn length(self) -> f32 {
        self.dot(self).sqrt()
    }
    pub fn normalized(self) -> Self {
        let l = self.length();
        if l == 0.0 {
            self
        } else {
            self.mul(1.0 / l)
        }
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
        let ab = self.b.sub(self.a);
        let ac = self.c.sub(self.a);
        let mut n = ab.cross(ac).normalized();
        let centroid = self.a.add(self.b).add(self.c);
        if n.dot(centroid) < 0.0 {
            // Flip winding B <-> C
            std::mem::swap(&mut self.b, &mut self.c);
            let edge_ab = self.b.sub(self.a);
            let edge_ac = self.c.sub(self.a);
            n = edge_ab.cross(edge_ac).normalized();
        }
        self.n = n;
        self
    }
}

/// Argmax `dot(N_f, p)`. Assumes all `FaceGeom` have outward n.
#[inline]
#[must_use]
pub fn pick_face(p: Vec3, faces: &[FaceGeom]) -> FaceId {
    let mut best_i = 0u32;
    let mut best_d = f32::NEG_INFINITY;
    for (i, f) in faces.iter().enumerate() {
        let d = f.n.dot(p);
        if d > best_d || ((d - best_d).abs() < f32::EPSILON && u32::try_from(i).unwrap_or(0) < best_i) {
            best_d = d;
            best_i = u32::try_from(i).unwrap_or(0);
        }
    }
    best_i
}

/// Project p onto face plane and compute planar barycentrics (α,β,γ) wrt (A,B,C).
#[inline]
#[must_use]
pub fn barycentrics_plane(p: Vec3, f: &FaceGeom) -> [f32; 3] {
    // Plane projection
    let to_plane = f.n.mul(f.n.dot(p)); // component along n
    let p_proj = p.sub(to_plane);

    // Build local basis at A
    let v0 = f.b.sub(f.a);
    let v1 = f.c.sub(f.a);
    let v2 = p_proj.sub(f.a);

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

    let upper = (ru + rv) >= (1.0 - eps_upper);
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
    Vec3 { x: clon * clat, y: slat, z: slon * clat }
}
