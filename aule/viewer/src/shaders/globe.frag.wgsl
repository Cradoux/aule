struct VSOut {
	@builtin(position) pos_clip : vec4<f32>,
	@location(0) world_pos : vec3<f32>,
	@location(1) height    : f32,
};
struct FSOut { @location(0) color : vec4<f32>; };

fn hypsometric(h: f32) -> vec3<f32> {
	let sea = vec3<f32>(0.15, 0.25, 0.75);
	let land = vec3<f32>(0.82, 0.78, 0.62);
	let t = clamp(h * 0.00015 + 0.5, 0.0, 1.0);
	return mix(sea, land, t);
}

@fragment
fn main(inf: VSOut) -> FSOut {
	var o : FSOut;
	let c = hypsometric(inf.height);
	o.color = vec4<f32>(c, 1.0);
	return o;
}


