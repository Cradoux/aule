struct Uniforms {
  width: u32,
  height: u32,
  F: u32,
  palette_mode: u32,
  debug_flags: u32,
  d_max: f32,
  h_max: f32,
  snowline: f32,
  eta_m: f32,
  inv_dmax: f32,
  inv_hmax: f32,
};

@group(0) @binding(0) var<uniform> U: Uniforms;
@group(0) @binding(1) var<storage, read> FACE_VERT_IDS: array<u32>;
@group(0) @binding(2) var<storage, read> FACE_OFFSETS: array<u32>;
@group(0) @binding(3) var lut_tex: texture_1d<f32>;
@group(0) @binding(4) var lut_smp: sampler;
@group(0) @binding(5) var<storage, read> FACE_GEOM: array<vec4<f32>>;
@group(0) @binding(6) var<storage, read> VERT_VALUE: array<f32>;
@group(0) @binding(7) var OUT_TEX: texture_storage_2d<rgba8unorm, write>;
// bindings 8 and 9 were previously triangle tables; no longer used
// binding(10): packed neighbor edge info: for each face (0..20), for edge idx 0..2 (opp A,B,C), store [neighbor_face, idx_shared0, idx_shared1]
@group(0) @binding(10) var<storage, read> FACE_EDGE_INFO: array<u32>;
// binding(11): per-pixel debug output: [face, tri_idx] for each pixel (x,y)
@group(0) @binding(11) var<storage, read_write> DEBUG_FACE_TRI: array<u32>;
// binding(12): optional CPU-provided per-pixel face pick (row-major), enabled by debug bit 7
@group(0) @binding(12) var<storage, read> CPU_FACE: array<u32>;
// binding(13): per-pixel raw (alpha,beta) written under debug bit 11
@group(0) @binding(13) var<storage, read_write> RAW_AB: array<vec2<f32>>;

const EPS_ROLLOVER: f32 = 1e-6;
const EPS_UPPER: f32 = 1e-6;
fn dot3(a: vec3<f32>, b: vec3<f32>) -> f32 { return a.x*b.x + a.y*b.y + a.z*b.z; }
fn cross3(a: vec3<f32>, b: vec3<f32>) -> vec3<f32> { return vec3<f32>(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
fn norm3(a: vec3<f32>) -> vec3<f32> { let l = max(length(a), 1e-9); return a / l; }
fn sph_area(a: vec3<f32>, b: vec3<f32>, c: vec3<f32>) -> f32 {
  let num = abs(dot3(cross3(a, b), c));
  let den = 1.0 + dot3(a, b) + dot3(b, c) + dot3(c, a);
  return 2.0 * atan2(num, max(den, 1e-9));
}

fn palette_color_from_lut(elev: f32) -> vec4<f32> {
  let e = elev;
  let idx = select(
    clamp((-e) * U.inv_dmax * 255.0, 0.0, 255.0),
    256.0 + clamp(e * U.inv_hmax * 255.0, 0.0, 255.0),
    e > 0.0
  );
  let u = (idx + 0.5) / 512.0;
  return textureSampleLevel(lut_tex, lut_smp, u, 0.0);
}

fn row_base(i: u32, F: u32) -> u32 { return i*(F + 1u) - (i*(i - 1u))/2u; }

fn face_pick(p: vec3<f32>) -> u32 {
  var best_f: u32 = 0u;
  var best_d: f32 = -1e9;
  for (var f: u32 = 0u; f < 20u; f = f + 1u) {
    let n = FACE_GEOM[4u*f + 3u].xyz;
    let d = dot3(n, p);
    if (d > best_d || (d == best_d && f < best_f)) { best_d = d; best_f = f; }
  }
  return best_f;
}

fn barycentric_in_face_raw(p: vec3<f32>, A: vec3<f32>, B: vec3<f32>, C: vec3<f32>) -> vec3<f32> {
  // Plane-projection barycentrics to match CPU exactly (raw, unclamped)
  let n = norm3(cross3(B - A, C - A));
  let t = dot3(n, p - A);
  let q = p - n * t;
  let v0 = B - A;
  let v1 = C - A;
  let v2 = q - A;
  let d00 = dot3(v0, v0);
  let d01 = dot3(v0, v1);
  let d11 = dot3(v1, v1);
  let d20 = dot3(v2, v0);
  let d21 = dot3(v2, v1);
  let denom = max(d00 * d11 - d01 * d01, 1e-12);
  let w_b = (d11 * d20 - d01 * d21) / denom;
  let w_c = (d00 * d21 - d01 * d20) / denom;
  let w_a = 1.0 - w_b - w_c;
  return vec3<f32>(w_a, w_b, w_c);
}

@compute @workgroup_size(8,8,1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  if (gid.x >= U.width || gid.y >= U.height) { return; }
  let fx = f32(gid.x) + 0.5;
  let fy = f32(gid.y) + 0.5;
  let lon = (fx / f32(U.width)) * 6.2831853 - 3.14159265;
  let lat = 1.57079633 - (fy / f32(U.height)) * 3.14159265;
  let cl = cos(lat);
  let p = vec3<f32>(cl * cos(lon), sin(lat), cl * sin(lon));

  var f: u32 = face_pick(p);
  if ((U.debug_flags & (1u<<7)) != 0u) {
    let lin = gid.y * U.width + gid.x;
    f = CPU_FACE[lin];
  }
  // If debug bit 10 set: override picked face with CPU-provided face, then recompute barycentrics
  if ((U.debug_flags & (1u<<10)) != 0u) {
    let lin = gid.y * U.width + gid.x;
    let f_cpu = CPU_FACE[lin];
    f = f_cpu;
  }
  // One-step neighbor rollover if outside
  var A = FACE_GEOM[4u*f + 0u].xyz;
  var B = FACE_GEOM[4u*f + 1u].xyz;
  var C = FACE_GEOM[4u*f + 2u].xyz;
  // Raw barycentrics for rollover decision
  var bc = barycentric_in_face_raw(p, A, B, C);
  let eps_roll = EPS_ROLLOVER;
  // One-hop neighbor rollover per spec
  for (var s: u32 = 0u; s < 1u; s = s + 1u) {
    if (bc.x >= -eps_roll && bc.y >= -eps_roll && bc.z >= -eps_roll) { break; }
    // pick edge: 0->opp A, 1->opp B, 2->opp C
    var edge: u32 = 0u;
    if (bc.x < bc.y && bc.x < bc.z) { edge = 0u; }
    else if (bc.y < bc.z) { edge = 1u; } else { edge = 2u; }
    let off = (f * 3u + edge) * 3u;
    let nf = FACE_EDGE_INFO[off + 0u];
    f = nf;
    A = FACE_GEOM[4u*f + 0u].xyz;
    B = FACE_GEOM[4u*f + 1u].xyz;
    C = FACE_GEOM[4u*f + 2u].xyz;
    // Recompute barycentrics in the neighbor face's basis (destination perspective)
    bc = barycentric_in_face_raw(p, A, B, C);
  }
  // Clamp tiny negatives and renormalize to keep inside
  let uu = max(bc.x, 0.0);
  let vv = max(bc.y, 0.0);
  let ww = max(bc.z, 0.0);
  let s_bc = max(uu + vv + ww, 1e-9);
  bc = vec3<f32>(uu/s_bc, vv/s_bc, ww/s_bc);
  // Write raw (alpha,beta) for CPU lattice when debug bit 11 is set
  if ((U.debug_flags & (1u<<11)) != 0u) {
    let lin = gid.y * U.width + gid.x;
    RAW_AB[lin] = vec2<f32>(bc.x, bc.y);
  }
  // Align with CPU lattice and face vertex table:
  // i (rows) follows α toward A (C→A), j (cols) follows β along AB.
  // Therefore u <- α*F (bc.x), v <- β*F (bc.y).
  let F = U.F;
  var u = clamp(bc.x * f32(F), 0.0, f32(F));
  var v = clamp(bc.y * f32(F), 0.0, f32(F));
  var iu: u32 = u32(floor(u));
  var iv: u32 = u32(floor(v));
  var fu = clamp(u - f32(iu), 1e-4, 1.0 - 1e-4);
  var fv = clamp(v - f32(iv), 1e-4, 1.0 - 1e-4);
  if (iu + iv > F - 1u) {
    let s = iu + iv - (F - 1u);
    if (fu > fv) { iu = iu - s; }
    else { iv = iv - s; }
    // Recompute fractional parts after reflection so weights are consistent
    fu = clamp(u - f32(iu), 1e-4, 1.0 - 1e-4);
    fv = clamp(v - f32(iv), 1e-4, 1.0 - 1e-4);
  }
  // Authoritative upper rule: diagonal belongs to lower; strict EPS (use raw residuals, not clamped)
  let upper = ((u - f32(iu)) + (v - f32(iv))) >= (1.0 - EPS_UPPER);
  // Analytic tri index per senior dev: tri_idx = face*F*F + v*(2F - v) + 2*u + (upper?1:0)
  let tri_local = u32(iv) * (2u * F - u32(iv)) + 2u * u32(iu) + select(0u, 1u, upper);
  let tri_idx   = f * F*F + tri_local;
  // Derive vertex ids directly from face vertex table
  let f_off = FACE_OFFSETS[f];
  let id_v0 = FACE_VERT_IDS[f_off + row_base(u32(iu), F) + u32(iv)];
  let id_v1 = FACE_VERT_IDS[f_off + row_base(u32(iu + 1), F) + u32(iv)];
  let id_v2 = FACE_VERT_IDS[f_off + row_base(u32(iu), F) + u32(iv + 1)];
  let id_v3 = FACE_VERT_IDS[f_off + row_base(u32(iu + 1), F) + u32(iv + 1)];
  let id0 = select(id_v0, id_v1, upper);
  let id1 = select(id_v1, id_v3, upper);
  let id2 = select(id_v2, id_v2, upper);
  var w0: f32;
  var w1: f32;
  var w2: f32;
  if (upper) {
    // Upper triangle: V1(i+1,j), V3(i+1,j+1), V2(i,j+1)
    // With u=alpha, v=beta → V1 gets fu, V2 gets fv
    w0 = fu;                 // V1
    w1 = 1.0 - fu - fv;      // V3
    w2 = fv;                 // V2
  } else {
    // Lower triangle: V0(i,j), V1(i+1,j), V2(i,j+1)
    // With u=alpha, v=beta → V1 gets fu, V2 gets fv
    w0 = 1.0 - fu - fv;      // V0
    w1 = fu;                 // V1
    w2 = fv;                 // V2
  }
  // Final clamp and renormalize for numerical robustness at seams
  {
    let ww0 = max(w0, 0.0);
    let ww1 = max(w1, 0.0);
    let ww2 = max(w2, 0.0);
    let s = max(ww0 + ww1 + ww2, 1e-9);
    w0 = ww0 / s; w1 = ww1 / s; w2 = ww2 / s;
  }
  let depth = w0 * VERT_VALUE[id0] + w1 * VERT_VALUE[id1] + w2 * VERT_VALUE[id2];
  let elev = U.eta_m - depth;
  var c = palette_color_from_lut(elev);
  // Debug wireframe: thin line where any bary weight near edge
  if ((U.debug_flags & 1u) != 0u) {
    let wire_w = 0.02;
    let bw = min(min(w0, w1), w2);
    if (bw < wire_w) { c = mix(c, vec4<f32>(0.0, 0.0, 0.0, 1.0), 0.6); }
  }
  // Face tint
  if ((U.debug_flags & 2u) != 0u) {
    let t = fract(sin(f32(f)*12.9898)*43758.5453);
    let tinted = mix(c.rgb, vec3<f32>(t, 1.0 - t, 0.5*t), 0.15);
    c = vec4<f32>(tinted, c.a);
  }
  // Grid overlay every 8 cells
  if ((U.debug_flags & 4u) != 0u) {
    let mod8u = (iu & 7u) == 0u;
    let mod8v = (iv & 7u) == 0u;
    if (mod8u || mod8v) { c = mix(c, vec4<f32>(1.0,1.0,1.0,1.0), 0.35); }
  }
  // Tri parity
  if ((U.debug_flags & 8u) != 0u) {
    if (upper) { let r = mix(c.rgb, vec3<f32>(1.0,0.0,0.0), 0.15); c = vec4<f32>(r, c.a); }
    else { let g = mix(c.rgb, vec3<f32>(0.0,1.0,0.0), 0.15); c = vec4<f32>(g, c.a); }
  }
  // Tri index stripes (bit 5)
  if ((U.debug_flags & (1u<<5)) != 0u) {
    let h = f32((tri_idx & 255u)) / 255.0; // 0..1
    let band = step(0.5, fract(h * 16.0));
    let tint = vec3<f32>(0.2 + 0.8*band, 0.2, 1.0 - 0.8*band);
    c = mix(c, vec4<f32>(tint, 1.0), 0.35);
  }
  // Optional: nearest-vertex sampling for diagnosis (bit 4)
  if ((U.debug_flags & 16u) != 0u) {
    // choose nearest of the tri's three vertices based on largest weight
    var wmax = w0; var vid = id0;
    if (w1 > wmax) { wmax = w1; vid = id1; }
    if (w2 > wmax) { wmax = w2; vid = id2; }
    let depth_nv = VERT_VALUE[vid];
    let elev_nv = U.eta_m - depth_nv;
    c = mix(c, palette_color_from_lut(elev_nv), 0.85);
  }
  // Face-id flat color (bit 9)
  if ((U.debug_flags & (1u<<9)) != 0u) {
    let h = f32(f & 255u) / 255.0;
    let c_flat = vec4<f32>(h, h * 0.5, 1.0 - h, 1.0);
    textureStore(OUT_TEX, vec2<i32>(i32(gid.x), i32(gid.y)), c_flat);
    return;
  }
  // Diagnostic write: face and tri_idx per pixel (guarded by bit 6)
  if ((U.debug_flags & (1u<<6)) != 0u) {
    let lin = gid.y * U.width + gid.x;
    let off = lin * 2u;
    DEBUG_FACE_TRI[off + 0u] = f;
    DEBUG_FACE_TRI[off + 1u] = tri_idx;
  }
  textureStore(OUT_TEX, vec2<i32>(i32(gid.x), i32(gid.y)), c);
}




