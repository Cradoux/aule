//! Exact pixel -> lon/lat mapping used by both CPU parity and the shader.

/// Map pixel center (x,y) and framebuffer size (width,height) to (lon, lat) in radians.
pub fn pixel_to_lon_lat(x: u32, y: u32, width: u32, height: u32) -> (f32, f32) {
	// pixel center in [0,1]
	let u = (x as f32 + 0.5) / (width as f32);
	let v = (y as f32 + 0.5) / (height as f32);
	// Equirectangular: lon in [-pi, pi], lat in [-pi/2, pi/2]
	// NOTE: v is top->bottom, so latitude uses (0.5 - v)
	let lon = core::f32::consts::TAU * (u - 0.5);
	let lat = core::f32::consts::PI * (0.5 - v);
	(lon, lat)
}

#[inline]
pub fn sph_to_unit(lon: f32, lat: f32) -> aule_geo::Vec3 {
	// Delegate to shared implementation to ensure parity with CPU/GPU
	aule_geo::sph_to_unit(lon, lat)
}
