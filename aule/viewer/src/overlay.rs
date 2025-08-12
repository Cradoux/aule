use egui::{epaint::Shape, Color32, Pos2, Rect, Stroke, Vec2};

pub struct OverlayState {
    pub show_hud: bool,
    pub show_plates: bool,
    pub show_vel: bool,
    pub show_bounds: bool,

    // HUD parameters
    pub vel_scale_px_per_cm_yr: f32, // px per (cm/yr)
    pub max_arrows_slider: u32,
    pub max_bounds_slider: u32,
    pub adaptive_cap: bool,
    pub live_arrows_cap: u32,
    pub live_bounds_cap: u32,

    // Cache management
    last_rect_key: Option<[f32; 4]>,
    last_vel_scale: f32,
    last_arrows_cap: u32,
    last_bounds_cap: u32,
    plates_cache: Option<Vec<Shape>>,
    vel_cache: Option<Vec<Shape>>,
    bounds_cache: Option<Vec<Shape>>,
}

impl Default for OverlayState {
    fn default() -> Self {
        let max_default: u32 = 4000;
        let vel_scale_default: f32 = 0.5; // px per (cm/yr)
        Self {
            show_hud: true,
            show_plates: false,
            show_vel: false,
            show_bounds: false,
            vel_scale_px_per_cm_yr: vel_scale_default,
            max_arrows_slider: max_default,
            max_bounds_slider: max_default,
            adaptive_cap: false,
            live_arrows_cap: max_default,
            live_bounds_cap: max_default,
            last_rect_key: None,
            last_vel_scale: vel_scale_default,
            last_arrows_cap: max_default,
            last_bounds_cap: max_default,
            plates_cache: None,
            vel_cache: None,
            bounds_cache: None,
        }
    }
}

pub fn color_for_plate(id: u16) -> Color32 {
    // Simple hash to color
    let h = (id as u32).wrapping_mul(0x9E37_79B9);
    let r = (h & 0xFF) as u8;
    let g = ((h >> 8) & 0xFF) as u8;
    let b = ((h >> 16) & 0xFF) as u8;
    Color32::from_rgb(r, g, b)
}

pub fn project_equirect(lat: f32, lon: f32, rect: Rect) -> Pos2 {
    let w = rect.width();
    let h = rect.height();
    let x = ((lon as f64 + std::f64::consts::PI) / (2.0 * std::f64::consts::PI)) as f32 * w;
    let y = (1.0 - ((lat as f64 + std::f64::consts::FRAC_PI_2) / std::f64::consts::PI) as f32) * h;
    Pos2::new(rect.left() + x, rect.top() + y)
}

pub fn build_plate_points(rect: Rect, latlon: &[[f32; 2]], plate_id: &[u16]) -> Vec<Shape> {
    let mut shapes = Vec::with_capacity(plate_id.len());
    for (i, &pid) in plate_id.iter().enumerate() {
        let p = project_equirect(latlon[i][0], latlon[i][1], rect);
        let col = color_for_plate(pid);
        shapes.push(Shape::circle_filled(p, 1.0, col));
    }
    shapes
}

pub fn build_velocity_arrows(
    rect: Rect,
    latlon: &[[f32; 2]],
    vel_en: &[[f32; 2]],
    max_arrows: usize,
    px_per_cm_yr: f32,
) -> Vec<Shape> {
    let n = latlon.len();
    let step = (n / max_arrows.max(1)).max(1);
    let mut shapes = Vec::new();
    for i in (0..n).step_by(step) {
        let p = project_equirect(latlon[i][0], latlon[i][1], rect);
        // Convert (m/yr) -> (cm/yr) and scale to pixels
        let v =
            Vec2::new(vel_en[i][0] * 100.0 * px_per_cm_yr, -vel_en[i][1] * 100.0 * px_per_cm_yr);
        shapes.push(Shape::line_segment([p, p + v], Stroke::new(1.0, Color32::WHITE)));
    }
    shapes
}

pub fn build_boundary_strokes(
    rect: Rect,
    latlon: &[[f32; 2]],
    edges: &[(u32, u32, u8)],
    max_edges: usize,
) -> Vec<Shape> {
    let step = (edges.len() / max_edges.max(1)).max(1);
    let mut shapes = Vec::new();
    for (idx, &(u, v, class)) in edges.iter().enumerate().step_by(step) {
        let pu = project_equirect(latlon[u as usize][0], latlon[u as usize][1], rect);
        let pv = project_equirect(latlon[v as usize][0], latlon[v as usize][1], rect);
        let center = Pos2::new((pu.x + pv.x) * 0.5, (pu.y + pv.y) * 0.5);
        let dir = (pv - pu).normalized();
        let orth = Vec2::new(-dir.y, dir.x);
        let half = orth * 3.0;
        let color = match class {
            1 => Color32::BLUE,
            2 => Color32::RED,
            3 => Color32::YELLOW,
            _ => Color32::GRAY,
        };
        shapes.push(Shape::line_segment([center - half, center + half], Stroke::new(1.5, color)));
        if idx >= edges.len() {
            break;
        }
    }
    shapes
}

fn rect_key(rect: Rect) -> [f32; 4] {
    [rect.left(), rect.top(), rect.right(), rect.bottom()]
}

impl OverlayState {
    pub fn effective_arrows_cap(&self) -> u32 {
        if self.adaptive_cap {
            self.live_arrows_cap
        } else {
            self.max_arrows_slider
        }
    }
    pub fn effective_bounds_cap(&self) -> u32 {
        if self.adaptive_cap {
            self.live_bounds_cap
        } else {
            self.max_bounds_slider
        }
    }

    pub fn ensure_params_and_invalidate_if_needed(
        &mut self,
        rect: Rect,
        arrows_cap: u32,
        bounds_cap: u32,
    ) {
        let rk = rect_key(rect);
        if self.last_rect_key.map_or(true, |k| k != rk) {
            self.last_rect_key = Some(rk);
            self.plates_cache = None;
            self.vel_cache = None;
            self.bounds_cache = None;
        }
        if (self.last_vel_scale - self.vel_scale_px_per_cm_yr).abs() > f32::EPSILON
            || self.last_arrows_cap != arrows_cap
        {
            self.last_vel_scale = self.vel_scale_px_per_cm_yr;
            self.last_arrows_cap = arrows_cap;
            self.vel_cache = None;
        }
        if self.last_bounds_cap != bounds_cap {
            self.last_bounds_cap = bounds_cap;
            self.bounds_cache = None;
        }
        // Clamp live caps into slider ranges in case sliders changed
        self.live_arrows_cap = self.live_arrows_cap.clamp(self.min_cap(), self.max_cap_arrows());
        self.live_bounds_cap = self.live_bounds_cap.clamp(self.min_cap(), self.max_cap_bounds());
    }

    fn min_cap(&self) -> u32 {
        500
    }
    fn max_cap_arrows(&self) -> u32 {
        self.max_arrows_slider.max(500)
    }
    fn max_cap_bounds(&self) -> u32 {
        self.max_bounds_slider.max(500)
    }

    pub fn shapes_for_plates(
        &mut self,
        rect: Rect,
        latlon: &[[f32; 2]],
        plate_id: &[u16],
    ) -> &[Shape] {
        if self.plates_cache.is_none() {
            self.plates_cache = Some(build_plate_points(rect, latlon, plate_id));
        }
        match self.plates_cache.as_deref() {
            Some(s) => s,
            None => &[],
        }
    }

    pub fn shapes_for_velocities(
        &mut self,
        rect: Rect,
        latlon: &[[f32; 2]],
        vel_en: &[[f32; 2]],
    ) -> &[Shape] {
        if self.vel_cache.is_none() {
            let cap = self.effective_arrows_cap() as usize;
            self.vel_cache =
                Some(build_velocity_arrows(rect, latlon, vel_en, cap, self.vel_scale_px_per_cm_yr));
        }
        match self.vel_cache.as_deref() {
            Some(s) => s,
            None => &[],
        }
    }

    pub fn shapes_for_boundaries(
        &mut self,
        rect: Rect,
        latlon: &[[f32; 2]],
        edges: &[(u32, u32, u8)],
    ) -> &[Shape] {
        if self.bounds_cache.is_none() {
            let cap = self.effective_bounds_cap() as usize;
            self.bounds_cache = Some(build_boundary_strokes(rect, latlon, edges, cap));
        }
        match self.bounds_cache.as_deref() {
            Some(s) => s,
            None => &[],
        }
    }

    pub fn update_adaptive_caps(&mut self, dt_ms: f32) {
        if !self.adaptive_cap {
            self.live_arrows_cap = self.max_arrows_slider.max(500);
            self.live_bounds_cap = self.max_bounds_slider.max(500);
            return;
        }
        // Simple ±10% adjustment per frame toward 16.6 ms target
        let target = 16.6f32;
        if dt_ms > target {
            let dec = |x: u32| ((x as f32) * 0.9).round() as u32;
            let new_ar = dec(self.live_arrows_cap).max(self.min_cap());
            let new_bd = dec(self.live_bounds_cap).max(self.min_cap());
            if new_ar != self.live_arrows_cap {
                self.vel_cache = None;
            }
            if new_bd != self.live_bounds_cap {
                self.bounds_cache = None;
            }
            self.live_arrows_cap = new_ar.min(self.max_cap_arrows());
            self.live_bounds_cap = new_bd.min(self.max_cap_bounds());
        } else {
            let inc = |x: u32| ((x as f32) * 1.1).round() as u32;
            let new_ar = inc(self.live_arrows_cap).min(self.max_cap_arrows());
            let new_bd = inc(self.live_bounds_cap).min(self.max_cap_bounds());
            if new_ar != self.live_arrows_cap {
                self.vel_cache = None;
            }
            if new_bd != self.live_bounds_cap {
                self.bounds_cache = None;
            }
            self.live_arrows_cap = new_ar.max(self.min_cap());
            self.live_bounds_cap = new_bd.max(self.min_cap());
        }
    }
}
