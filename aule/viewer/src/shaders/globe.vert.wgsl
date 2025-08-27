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
@group(0) @binding(0) var<uniform> G : Globals;
@group(1) @binding(0) var<storage, read> ZBuf : array<f32>;

struct VSIn {
	@location(0) pos_unit : vec3<f32>,
	@location(1) vid      : u32,
};
struct VSOut {
	@builtin(position) pos_clip : vec4<f32>,
	@location(0) world_pos : vec3<f32>,
	@location(1) height    : f32,
};

@vertex
fn main(input: VSIn) -> VSOut {
	let u = normalize(input.pos_unit);
	let h = ZBuf[input.vid];
	// Clamp radius to avoid inward flips when exaggeration is large negative or huge positive
	let r = G.radius + clamp(h * G.exagger, -0.2, 0.2);
	let world = u * r;
	var o : VSOut;
	o.pos_clip = G.view_proj * vec4<f32>(world, 1.0);
	o.world_pos = world;
	o.height = h;
	return o;
}


