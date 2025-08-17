struct Uniforms {
  width: u32,
  height: u32,
  F: u32,
  palette_mode: u32,
  _pad0: u32,
  d_max: f32,
  h_max: f32,
  snowline: f32,
  eta_m: f32,
};

@group(0) @binding(0) var<uniform> U: Uniforms;
@group(0) @binding(1) var<storage, read> FACE_VERT_IDS: array<u32>;
@group(0) @binding(2) var<storage, read> FACE_OFFSETS: array<u32>;
@group(0) @binding(3) var<storage, read> FACE_GEOM: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read> VERT_VALUE: array<f32>;
@group(0) @binding(5) var OUT_TEX: texture_storage_2d<rgba8unorm, write>;

fn dot3(a: vec3<f32>, b: vec3<f32>) -> f32 { return a.x*b.x + a.y*b.y + a.z*b.z; }
fn cross3(a: vec3<f32>, b: vec3<f32>) -> vec3<f32> { return vec3<f32>(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
fn norm3(a: vec3<f32>) -> vec3<f32> { let l = max(length(a), 1e-9); return a / l; }

fn color_hypso(val: f32, dmax: f32, hmax: f32, snow: f32) -> vec4<f32> {
  if (val >= 0.0) {
    let t = clamp(val / max(dmax, 1.0), 0.0, 1.0);
    let c0 = vec3<f32>(0.024, 0.137, 0.259);
    let c1 = vec3<f32>(0.110, 0.400, 0.639);
    let c2 = vec3<f32>(0.459, 0.737, 0.824);
    let c = mix(c0, c1, min(t, 0.7)/0.7);
    let c3 = mix(c, c2, max(t-0.7, 0.0)/0.3);
    return vec4<f32>(c3, 1.0);
  } else {
    let h = -val;
    if (h >= snow) {
      let t = clamp((h - snow) / max(hmax - snow, 1.0), 0.0, 1.0);
      let high = vec3<f32>(0.73, 0.75, 0.78);
      let snowc= vec3<f32>(0.96, 0.98, 1.00);
      return vec4<f32>(mix(high, snowc, t), 1.0);
    }
    let t = clamp(h / max(snow, 1.0), 0.0, 1.0);
    let beach = vec3<f32>(0.925, 0.878, 0.729);
    let hill  = vec3<f32>(0.573, 0.506, 0.392);
    let high  = vec3<f32>(0.73,  0.75,  0.78);
    let c = select(mix(beach, hill, t/0.5), mix(hill, high, (t-0.5)/0.5), t > 0.5);
    return vec4<f32>(c, 1.0);
  }
}

fn palette_color(val: f32) -> vec4<f32> {
  if (U.palette_mode == 0u) {
    return color_hypso(val, U.d_max, U.h_max, U.snowline);
  } else {
    return color_hypso(val, U.d_max, U.h_max, U.snowline);
  }
}

fn row_base(i: u32, F: u32) -> u32 { return i*(F + 1u) - (i*(i - 1u))/2u; }

fn face_pick(p: vec3<f32>) -> u32 {
  var best_f: u32 = 0u;
  var best_d: f32 = -1e9;
  for (var f: u32 = 0u; f < 20u; f = f + 1u) {
    let n = FACE_GEOM[4u*f + 3u].xyz;
    let d = dot3(n, p);
    if (d > best_d) { best_d = d; best_f = f; }
  }
  return best_f;
}

fn barycentric_in_face(p: vec3<f32>, A: vec3<f32>, B: vec3<f32>, C: vec3<f32>) -> vec3<f32> {
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
  let v = (d11 * d20 - d01 * d21) / denom;
  let w = (d00 * d21 - d01 * d20) / denom;
  let u = 1.0 - v - w;
  let uu = max(u, 0.0);
  let vv = max(v, 0.0);
  let ww = max(w, 0.0);
  let s = max(uu + vv + ww, 1e-9);
  return vec3<f32>(uu/s, vv/s, ww/s);
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

  let f = face_pick(p);
  let A = FACE_GEOM[4u*f + 0u].xyz;
  let B = FACE_GEOM[4u*f + 1u].xyz;
  let C = FACE_GEOM[4u*f + 2u].xyz;
  let bc = barycentric_in_face(p, A, B, C);
  let F = U.F;
  var u = clamp(bc.y * f32(F), 0.0, f32(F));
  var v = clamp(bc.z * f32(F), 0.0, f32(F));
  var i = u32(floor(u));
  var j = u32(floor(v));
  if (i + j > F) {
    let excess = i + j - F;
    if (j >= excess) { j -= excess; } else { i -= excess - j; j = 0u; }
  }
  var fu = clamp(u - f32(i), 1e-4, 1.0 - 1e-4);
  var fv = clamp(v - f32(j), 1e-4, 1.0 - 1e-4);

  let base = FACE_OFFSETS[f] + row_base(i, F) + j;
  var id0: u32;
  var id1: u32;
  var id2: u32;
  var w0: f32;
  var w1: f32;
  var w2: f32;
  if (fu + fv <= 1.0) {
    let base_i1 = FACE_OFFSETS[f] + row_base(i+1u, F);
    id0 = FACE_VERT_IDS[ base ];
    id1 = FACE_VERT_IDS[ base_i1 + j ];
    id2 = FACE_VERT_IDS[ base + 1u ];
    w0 = 1.0 - fu - fv; w1 = fu; w2 = fv;
  } else {
    let base_i1 = FACE_OFFSETS[f] + row_base(i+1u, F);
    id0 = FACE_VERT_IDS[ base_i1 + j ];
    id1 = FACE_VERT_IDS[ base_i1 + (j + 1u) ];
    id2 = FACE_VERT_IDS[ base + 1u ];
    let fu2 = 1.0 - fu; let fv2 = 1.0 - fv;
    w0 = 1.0 - fu2 - fv2; w1 = fv2; w2 = fu2;
  }

  let depth = w0 * VERT_VALUE[id0] + w1 * VERT_VALUE[id1] + w2 * VERT_VALUE[id2];
  let elev = U.eta_m - depth;
  let c = palette_color(elev);
  textureStore(OUT_TEX, vec2<i32>(i32(gid.x), i32(gid.y)), c);
}




