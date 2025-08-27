struct VSOut {
	@builtin(position) pos_clip : vec4<f32>,
	@location(0) world_pos : vec3<f32>,
	@location(1) height    : f32,
};
struct FSOut { @location(0) color : vec4<f32> };

@group(0) @binding(1) var lut_tex : texture_1d<f32>;
@group(0) @binding(2) var lut_samp : sampler;
@group(0) @binding(0) var<uniform> G : Globals;

struct Globals {
    view_proj : mat4x4<f32>,
    radius    : f32,
    exagger   : f32,
    debug_flags : u32,
    _pad : f32,
    d_max : f32,
    h_max : f32,
    _pad2 : vec2<f32>,
};

@fragment
fn main(inf: VSOut) -> FSOut {
	var o : FSOut;
	let z = inf.height;
	var idx : f32;
	// Split at 0 (sea level already applied to z)
	if (z >= 0.0) {
		let t = clamp(z / max(1e-6, G.h_max), 0.0, 1.0);
		idx = 256.0 + t * 255.0;
	} else {
		let t = clamp((-z) / max(1e-6, G.d_max), 0.0, 1.0);
		idx = t * 255.0;
	}
	let u = (idx + 0.5) / 512.0;
	let c = textureSample(lut_tex, lut_samp, u).rgb;
	o.color = vec4<f32>(c, 1.0);
	return o;
}


