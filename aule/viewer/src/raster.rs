use egui::ColorImage;

use crate::overlay::OverlayState;
use engine::world::World;

pub fn render_map(world: &World, ov: &OverlayState, width: u32, height: u32) -> ColorImage {
    let w = width as usize;
    let h = height as usize;
    let mut img = ColorImage::new([w, h], egui::Color32::BLACK);

    // Per-vertex scalar (depth per grid vertex)
    let values: &[f32] = &world.depth_m;

    // Precompute sea-aware scales
    let d_max = ov.hypso_d_max.max(1.0);
    let h_max = ov.hypso_h_max.max(1.0);
    let snow = ov.hypso_snowline;

    let t_start = std::time::Instant::now();

    // Build a coarse lat/lon bin index to avoid O(Ncells) nearest search
    let bins_x = (width as usize).max(256);
    let bins_y = (height as usize / 2).max(128);
    let mut bin: Vec<i32> = vec![-1; bins_x * bins_y];
    for (i, ll) in world.grid.latlon.iter().enumerate() {
        let lat = ll[0];
        let lon = ll[1];
        let u = ((lon + std::f32::consts::PI) / std::f32::consts::TAU).clamp(0.0, 0.9999999);
        let v = (0.5 - lat / std::f32::consts::PI).clamp(0.0, 0.9999999);
        let x = (u * (bins_x as f32)) as usize;
        let y = (v * (bins_y as f32)) as usize;
        bin[y * bins_x + x] = i as i32;
    }

    for y in 0..h {
        let v = (y as f32 + 0.5) / (h as f32);
        let lat = std::f32::consts::FRAC_PI_2 - v * std::f32::consts::PI;
        let sin_lat = lat.sin();
        let cos_lat = lat.cos();
        for x in 0..w {
            let u = (x as f32 + 0.5) / (w as f32);
            let lon = u * std::f32::consts::TAU - std::f32::consts::PI;
            let sin_lon = lon.sin();
            let cos_lon = lon.cos();
            let _p = [cos_lat * cos_lon, sin_lat, cos_lat * sin_lon];

            // Coarse bin lookup with 3x3 neighborhood search
            let bx =
                (((lon + std::f32::consts::PI) / std::f32::consts::TAU) * (bins_x as f32)) as isize;
            let by = ((0.5 - lat / std::f32::consts::PI) * (bins_y as f32)) as isize;
            let mut best_i = 0usize;
            let mut best_d = f32::INFINITY;
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let xx = (bx + dx).rem_euclid(bins_x as isize) as usize;
                    let yy = ((by + dy).clamp(0, (bins_y as isize) - 1)) as usize;
                    let idx = bin[yy * bins_x + xx];
                    if idx >= 0 {
                        let i = idx as usize;
                        let ll = world.grid.latlon[i];
                        let dlat = lat - ll[0];
                        let mut dlon = lon - ll[1];
                        if dlon > std::f32::consts::PI {
                            dlon -= std::f32::consts::TAU;
                        }
                        if dlon < -std::f32::consts::PI {
                            dlon += std::f32::consts::TAU;
                        }
                        let dsq = dlat * dlat + dlon * dlon * cos_lat * cos_lat;
                        if dsq < best_d {
                            best_d = dsq;
                            best_i = i;
                        }
                    }
                }
            }
            let val = values[best_i];
            let col = if val >= 0.0 {
                // Ocean
                super::overlay::ocean_color32(val, d_max)
            } else {
                super::overlay::land_color32(-val, h_max, snow)
            };
            img.pixels[y * w + x] = col;
        }
    }

    let dt = t_start.elapsed().as_millis();
    println!(
        "[viewer] raster build W={} H={} | F={} | time={} ms",
        width, height, world.grid.frequency, dt
    );
    img
}
