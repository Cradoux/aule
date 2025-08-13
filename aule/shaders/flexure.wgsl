// Thin-plate + Winkler operator A w = D (D4x + D4y) + k w
// Memory layout: row-major, including halos; params carry total dims

struct Params {
    w_tot: u32,
    h_tot: u32,
    halo: u32,
    _pad: u32,
    dx: f32,
    d: f32,
    k: f32,
    _pad2: f32,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> W: array<f32>;
@group(0) @binding(2) var<storage, read_write> OUT: array<f32>;

fn idx(i: u32, j: u32) -> u32 { return j * P.w_tot + i; }

@compute @workgroup_size(8, 8, 1)
fn apply_A(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x + P.halo;
    let j = gid.y + P.halo;
    if (i >= P.w_tot - P.halo || j >= P.h_tot - P.halo) { return; }
    let inv_dx4 = 1.0 / (P.dx * P.dx * P.dx * P.dx);
    let c = idx(i, j);
    let w00 = W[c];
    // x 1D 4th derivative
    let wxm2 = W[idx(i - 2u, j)];
    let wxm1 = W[idx(i - 1u, j)];
    let wxp1 = W[idx(i + 1u, j)];
    let wxp2 = W[idx(i + 2u, j)];
    let d4x = (wxm2 - 4.0 * wxm1 + 6.0 * w00 - 4.0 * wxp1 + wxp2) * inv_dx4;
    // y 1D 4th derivative
    let wym2 = W[idx(i, j - 2u)];
    let wym1 = W[idx(i, j - 1u)];
    let wyp1 = W[idx(i, j + 1u)];
    let wyp2 = W[idx(i, j + 2u)];
    let d4y = (wym2 - 4.0 * wym1 + 6.0 * w00 - 4.0 * wyp1 + wyp2) * inv_dx4;
    OUT[c] = P.d * (d4x + d4y) + P.k * w00;
}


