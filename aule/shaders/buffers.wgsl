// T-020: WGSL buffer declarations and trivial identity compute for tests.

struct F32Buf { data: array<f32>; };
struct U32Buf { data: array<u32>; };

@group(0) @binding(0) var<storage, read_write> h: F32Buf;
@group(0) @binding(1) var<storage, read_write> th_c: F32Buf;
@group(0) @binding(2) var<storage, read_write> age_ocean: F32Buf;
@group(0) @binding(3) var<storage, read_write> S: F32Buf;
@group(0) @binding(4) var<storage, read_write> V: F32Buf;
@group(0) @binding(5) var<storage, read_write> plate_id: U32Buf;
@group(0) @binding(6) var<storage, read_write> B: F32Buf;

@compute @workgroup_size(64)
fn cs_identity(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    // trivial touches to prove bindings; identity op
    h.data[i] = h.data[i];
    th_c.data[i] = th_c.data[i];
    age_ocean.data[i] = age_ocean.data[i];
    S.data[i] = S.data[i];
    V.data[i] = V.data[i];
    plate_id.data[i] = plate_id.data[i];
    B.data[i] = B.data[i];
}


