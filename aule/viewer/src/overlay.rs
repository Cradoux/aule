use egui::{
    epaint::{Mesh, Shape, Vertex, WHITE_UV},
    Color32, Pos2, Rect, Stroke, Vec2,
};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct ContKey {
    pub seed: u64,
    pub n: u32,
    pub radius_km: u32,
    pub falloff_km: u32,
    pub f: u32,
}

#[allow(dead_code)]
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
    last_subd_cap: u32,
    last_trans_cap: u32,
    pub(crate) plates_cache: Option<Vec<Shape>>,
    pub(crate) vel_cache: Option<Vec<Shape>>,
    pub(crate) bounds_cache: Option<Vec<Shape>>,

    // Age/bathymetry layers
    pub show_age: bool,
    pub show_bathy: bool,
    pub v_floor_cm_per_yr: f32,
    pub age_minmax: (f32, f32),
    pub depth_minmax: (f32, f32),
    pub lock_bathy_scale: bool,
    pub bathy_min_max: (f32, f32),
    pub target_ocean_fraction: f32,
    pub extra_offset_m: f32,
    pub apply_sea_level: bool,
    pub last_isostasy_offset_m: f64,
    pub(crate) age_cache: Option<Vec<Shape>>,
    pub(crate) bathy_cache: Option<Vec<Shape>>,

    // Age–Depth plot controls
    pub show_age_depth: bool,
    pub plot_sample_cap: u32,
    pub plot_bin_width_myr: f32,

    // Subduction bands
    pub show_subduction: bool,
    pub subd_trench: Option<Vec<Shape>>,
    pub subd_arc: Option<Vec<Shape>>,
    pub subd_backarc: Option<Vec<Shape>>,
    pub max_subd_slider: u32,
    pub live_subd_cap: u32,

    // Subduction params (viewer-controlled)
    pub _rollback_offset_km: f32,
    pub _rollback_rate_km_per_myr: f32,
    pub _backarc_extension_mode: bool,
    pub _backarc_extension_deepen_m: f32,

    // Transforms
    pub show_transforms: bool,
    pub trans_pull: Option<Vec<Shape>>,
    pub trans_rest: Option<Vec<Shape>>,
    pub trans_tau_open_m_per_yr: f32,
    pub trans_min_tangential_m_per_yr: f32,
    pub trans_basin_half_width_km: f32,
    pub trans_basin_deepen_m: f32,
    pub trans_ridge_like_uplift_m: f32,
    pub trans_max_points: u32,
    pub trans_pull_count: u32,
    pub trans_rest_count: u32,

    // Continents overlay/state
    pub show_continents: bool,
    pub cont_seed: u64,
    pub cont_n: u32,
    pub cont_radius_km: f64,
    pub cont_falloff_km: f64,
    pub cont_auto_amp: bool,
    pub cont_target_land_frac: f32,
    pub cont_manual_amp_m: f32,
    pub cont_max_points: usize,
    pub mesh_continents: Option<Mesh>,
    pub mesh_coastline: Option<Mesh>,
    pub cont_land_frac: f64,
    pub cont_amp_applied_m: f32,
    pub cont_key: Option<ContKey>,
    pub cont_template: Option<Vec<f32>>,

    // Flexure overlay/state
    pub show_flexure: bool,   // draw w overlay
    pub enable_flexure: bool, // apply w to depth
    pub E_gpa: f32,
    pub nu: f32,
    pub Te_km: f32,
    pub k_winkler: f32, // N/m^3
    pub wj_omega: f32,  // 0.6..0.9
    pub nu1: u32,
    pub nu2: u32,
    pub levels: u32,
    pub max_points_flex: usize, // overlay cap
    pub last_residual: f32,     // r_out / max(1,r_in)
    pub flex_mesh: Option<Mesh>,
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
            last_subd_cap: 10_000,
            last_trans_cap: 10_000,
            plates_cache: None,
            vel_cache: None,
            bounds_cache: None,
            show_age: false,
            show_bathy: false,
            v_floor_cm_per_yr: 0.5,
            age_minmax: (0.0, 0.0),
            depth_minmax: (0.0, 0.0),
            lock_bathy_scale: false,
            bathy_min_max: (2600.0, 6000.0),
            target_ocean_fraction: 0.70,
            extra_offset_m: 0.0,
            apply_sea_level: false,
            last_isostasy_offset_m: 0.0,
            age_cache: None,
            bathy_cache: None,
            show_age_depth: false,
            plot_sample_cap: 5000,
            plot_bin_width_myr: 5.0,
            show_subduction: false,
            subd_trench: None,
            subd_arc: None,
            subd_backarc: None,
            max_subd_slider: 10_000,
            live_subd_cap: 10_000,
            _rollback_offset_km: 0.0,
            _rollback_rate_km_per_myr: 0.0,
            _backarc_extension_mode: false,
            _backarc_extension_deepen_m: 600.0,
            show_transforms: false,
            trans_pull: None,
            trans_rest: None,
            trans_tau_open_m_per_yr: 0.01,
            trans_min_tangential_m_per_yr: 0.002,
            trans_basin_half_width_km: 60.0,
            trans_basin_deepen_m: 400.0,
            trans_ridge_like_uplift_m: -300.0,
            trans_max_points: 10_000,
            trans_pull_count: 0,
            trans_rest_count: 0,

            show_continents: false,
            cont_seed: 1_234_567,
            cont_n: 3,
            cont_radius_km: 2200.0,
            cont_falloff_km: 600.0,
            cont_auto_amp: true,
            cont_target_land_frac: 0.29,
            cont_manual_amp_m: 2500.0,
            cont_max_points: 20_000,
            mesh_continents: None,
            mesh_coastline: None,
            cont_land_frac: 0.0,
            cont_amp_applied_m: 0.0,
            cont_key: None,
            cont_template: None,

            show_flexure: false,
            enable_flexure: false,
            E_gpa: 70.0,
            nu: 0.25,
            Te_km: 25.0,
            k_winkler: 3.0e8f32,
            wj_omega: 0.7,
            nu1: 1,
            nu2: 1,
            levels: 2,
            max_points_flex: 10_000,
            last_residual: 0.0,
            flex_mesh: None,
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
    pub fn effective_subd_cap(&self) -> u32 {
        if self.adaptive_cap {
            self.live_subd_cap
        } else {
            self.max_subd_slider
        }
    }

    pub fn ensure_params_and_invalidate_if_needed(
        &mut self,
        rect: Rect,
        arrows_cap: u32,
        bounds_cap: u32,
        subd_cap: u32,
    ) {
        let rk = rect_key(rect);
        if self.last_rect_key.map_or(true, |k| k != rk) {
            self.last_rect_key = Some(rk);
            self.plates_cache = None;
            self.vel_cache = None;
            self.bounds_cache = None;
            self.subd_trench = None;
            self.subd_arc = None;
            self.subd_backarc = None;
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
        if self.last_subd_cap != subd_cap {
            self.last_subd_cap = subd_cap;
            self.subd_trench = None;
            self.subd_arc = None;
            self.subd_backarc = None;
            self.trans_pull = None;
            self.trans_rest = None;
        }
        if self.last_trans_cap != self.trans_max_points {
            self.last_trans_cap = self.trans_max_points;
            self.trans_pull = None;
            self.trans_rest = None;
        }
        // Clamp live caps into slider ranges in case sliders changed
        self.live_arrows_cap = self.live_arrows_cap.clamp(self.min_cap(), self.max_cap_arrows());
        self.live_bounds_cap = self.live_bounds_cap.clamp(self.min_cap(), self.max_cap_bounds());
        self.live_subd_cap = self.live_subd_cap.clamp(self.min_cap(), self.max_cap_subd());
        // Age/bathy caches are invalidated by rect changes (recompute positions)
        if self.last_rect_key.is_none() {
            self.age_cache = None;
            self.bathy_cache = None;
        }
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
    fn max_cap_subd(&self) -> u32 {
        self.max_subd_slider.max(500)
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
        self.plates_cache.as_deref().unwrap_or_default()
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
        self.vel_cache.as_deref().unwrap_or_default()
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
        self.bounds_cache.as_deref().unwrap_or_default()
    }

    pub fn update_adaptive_caps(&mut self, dt_ms: f32) {
        if !self.adaptive_cap {
            self.live_arrows_cap = self.max_arrows_slider.max(500);
            self.live_bounds_cap = self.max_bounds_slider.max(500);
            self.live_subd_cap = self.max_subd_slider.max(500);
            return;
        }
        // Simple ±10% adjustment per frame toward 16.6 ms target
        let target = 16.6f32;
        if dt_ms > target {
            let dec = |x: u32| ((x as f32) * 0.9).round() as u32;
            let new_ar = dec(self.live_arrows_cap).max(self.min_cap());
            let new_bd = dec(self.live_bounds_cap).max(self.min_cap());
            let new_sd = dec(self.live_subd_cap).max(self.min_cap());
            if new_ar != self.live_arrows_cap {
                self.vel_cache = None;
            }
            if new_bd != self.live_bounds_cap {
                self.bounds_cache = None;
            }
            if new_sd != self.live_subd_cap {
                self.subd_trench = None;
                self.subd_arc = None;
                self.subd_backarc = None;
            }
            self.live_arrows_cap = new_ar.min(self.max_cap_arrows());
            self.live_bounds_cap = new_bd.min(self.max_cap_bounds());
            self.live_subd_cap = new_sd.min(self.max_cap_subd());
        } else {
            let inc = |x: u32| ((x as f32) * 1.1).round() as u32;
            let new_ar = inc(self.live_arrows_cap).min(self.max_cap_arrows());
            let new_bd = inc(self.live_bounds_cap).min(self.max_cap_bounds());
            let new_sd = inc(self.live_subd_cap).min(self.max_cap_subd());
            if new_ar != self.live_arrows_cap {
                self.vel_cache = None;
            }
            if new_bd != self.live_bounds_cap {
                self.bounds_cache = None;
            }
            if new_sd != self.live_subd_cap {
                self.subd_trench = None;
                self.subd_arc = None;
                self.subd_backarc = None;
            }
            self.live_arrows_cap = new_ar.max(self.min_cap());
            self.live_bounds_cap = new_bd.max(self.min_cap());
            self.live_subd_cap = new_sd.max(self.min_cap());
        }
    }
}

#[inline]
fn sat01(x: f32) -> f32 {
    x.clamp(0.0, 1.0)
}

/// Viridis-like ramp
pub fn viridis_color32(t: f32) -> Color32 {
    const C: [(u8, u8, u8); 9] = [
        (68, 1, 84),
        (71, 44, 122),
        (59, 81, 139),
        (44, 113, 142),
        (33, 144, 141),
        (39, 173, 129),
        (92, 200, 99),
        (170, 220, 50),
        (253, 231, 37),
    ];
    let t = sat01(t);
    let segs = (C.len() - 1) as f32;
    let x = t * segs;
    let i = x.floor() as usize;
    let f = x - (i as f32);
    if i >= C.len() - 1 {
        let (r, g, b) = C[C.len() - 1];
        return Color32::from_rgb(r, g, b);
    }
    let (r0, g0, b0) = C[i];
    let (r1, g1, b1) = C[i + 1];
    let r = (r0 as f32 + f * ((r1 as f32) - (r0 as f32))).round() as u8;
    let g = (g0 as f32 + f * ((g1 as f32) - (g0 as f32))).round() as u8;
    let b = (b0 as f32 + f * ((b1 as f32) - (b0 as f32))).round() as u8;
    Color32::from_rgb(r, g, b)
}

pub fn viridis_map(value: f32, min_v: f32, max_v: f32) -> Color32 {
    let t = if max_v > min_v { ((value - min_v) / (max_v - min_v)).clamp(0.0, 1.0) } else { 0.0 };
    viridis_color32(t)
}

fn blue_ramp(value: f32, min_v: f32, max_v: f32) -> Color32 {
    let t = if max_v > min_v { ((value - min_v) / (max_v - min_v)).clamp(0.0, 1.0) } else { 0.0 };
    let r = (10.0 + 20.0 * t) as u8;
    let g = (30.0 + 80.0 * t) as u8;
    let b = (80.0 + 160.0 * t) as u8;
    Color32::from_rgb(r, g, b)
}

impl OverlayState {
    pub fn rebuild_age_shapes(&mut self, rect: Rect, latlon: &[[f32; 2]], age_myr: &[f32]) {
        let mut shapes = Vec::with_capacity(latlon.len());
        let (amin, amax) = self.age_minmax;
        for i in 0..latlon.len() {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            let col = viridis_map(age_myr[i], amin, amax);
            shapes.push(Shape::circle_filled(p, 1.2, col));
        }
        self.age_cache = Some(shapes);
    }

    pub fn rebuild_bathy_shapes(&mut self, rect: Rect, grid: &engine::grid::Grid, depth_m: &[f32]) {
        let mut shapes = Vec::with_capacity(grid.latlon.len());
        let (dmin, dmax) = self.depth_minmax;
        for (i, &ll) in grid.latlon.iter().enumerate() {
            let p = project_equirect(ll[0], ll[1], rect);
            let col = blue_ramp(depth_m[i], dmin, dmax);
            shapes.push(Shape::circle_filled(p, 1.2, col));
        }
        // Draw 0 m coastline as thin black segments where any neighbor crosses zero
        let mut coast = Vec::new();
        for i in 0..grid.latlon.len() {
            let di = depth_m[i];
            for &n in &grid.n1[i] {
                let j = n as usize;
                let dj = depth_m[j];
                if (di > 0.0 && dj <= 0.0) || (di <= 0.0 && dj > 0.0) {
                    let pu = project_equirect(grid.latlon[i][0], grid.latlon[i][1], rect);
                    let pv = project_equirect(grid.latlon[j][0], grid.latlon[j][1], rect);
                    let center = Pos2::new((pu.x + pv.x) * 0.5, (pu.y + pv.y) * 0.5);
                    let dir = (pv - pu).normalized();
                    let orth = Vec2::new(-dir.y, dir.x);
                    let half = orth * 1.0;
                    coast.push(Shape::line_segment(
                        [center - half, center + half],
                        Stroke::new(1.0, Color32::BLACK),
                    ));
                }
            }
        }
        shapes.extend(coast);
        self.bathy_cache = Some(shapes);
    }

    /// Build flexure deflection overlay as a lightweight point mesh colored by value.
    pub fn rebuild_flexure_mesh(
        &mut self,
        rect: Rect,
        grid: &engine::grid::Grid,
        w_m: &[f32],
        cap_total: usize,
    ) {
        let mut mesh = Mesh::default();
        let n = grid.latlon.len().min(w_m.len());
        let mut idxs: Vec<usize> = (0..n).collect();
        // Subsample uniformly
        let stride = (n / cap_total.max(1)).max(1);
        // Scale by min/max of w for color
        let mut wmin = f32::INFINITY;
        let mut wmax = f32::NEG_INFINITY;
        for &v in w_m.iter().take(n) {
            if v.is_finite() {
                if v < wmin {
                    wmin = v;
                }
                if v > wmax {
                    wmax = v;
                }
            }
        }
        if !wmin.is_finite() || !wmax.is_finite() || wmin >= wmax {
            wmin = 0.0;
            wmax = 1.0;
        }
        for &i in idxs.iter().step_by(stride) {
            let ll = grid.latlon[i];
            let p = project_equirect(ll[0], ll[1], rect);
            let col = viridis_map(w_m[i], wmin, wmax);
            Self::mesh_add_dot(&mut mesh, p, 1.2, col);
        }
        if mesh.vertices.is_empty() {
            self.flex_mesh = None;
        } else {
            self.flex_mesh = Some(mesh);
        }
    }

    /// Build continent land mask and coastline meshes
    pub fn rebuild_continent_meshes(
        &mut self,
        rect: Rect,
        grid: &engine::grid::Grid,
        depth_m: &[f32],
        mask_land: &[bool],
        cap_total: usize,
    ) -> (Mesh, Mesh) {
        // Land point cloud
        let mut land_indices: Vec<usize> = Vec::new();
        for (i, &is_land) in mask_land.iter().enumerate() {
            if is_land {
                land_indices.push(i);
            }
        }
        let mut mesh_land = Mesh::default();
        let mut mesh_coast = Mesh::default();
        let total = land_indices.len().max(1);
        let stride = (total / cap_total.max(1)).max(1);
        let col_land = Color32::from_rgba_unmultiplied(120, 120, 120, 160);
        let r = 1.2;
        for &i in land_indices.iter().step_by(stride) {
            let p = project_equirect(grid.latlon[i][0], grid.latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_land, p, r, col_land);
        }
        // Coastline: neighbors crossing 0 m
        for i in 0..grid.latlon.len() {
            let di = depth_m[i];
            for &n in &grid.n1[i] {
                let j = n as usize;
                let dj = depth_m[j];
                if (di > 0.0 && dj <= 0.0) || (di <= 0.0 && dj > 0.0) {
                    let pu = project_equirect(grid.latlon[i][0], grid.latlon[i][1], rect);
                    let pv = project_equirect(grid.latlon[j][0], grid.latlon[j][1], rect);
                    let center = Pos2::new((pu.x + pv.x) * 0.5, (pu.y + pv.y) * 0.5);
                    let dir = (pv - pu).normalized();
                    let _orth = Vec2::new(-dir.y, dir.x);
                    // Approximate thin segment by a small dot at the midpoint for performance
                    Self::mesh_add_dot(&mut mesh_coast, center, 0.8, Color32::BLACK);
                }
            }
        }
        (mesh_land, mesh_coast)
    }

    pub fn age_shapes(&self) -> &[Shape] {
        self.age_cache.as_deref().unwrap_or_default()
    }
    pub fn bathy_shapes(&self) -> &[Shape] {
        self.bathy_cache.as_deref().unwrap_or_default()
    }

    fn mesh_add_dot(mesh: &mut Mesh, p: Pos2, r: f32, color: Color32) {
        let base = mesh.vertices.len() as u32;
        let v0 = Vertex { pos: [p.x - r, p.y - r].into(), uv: WHITE_UV, color };
        let v1 = Vertex { pos: [p.x + r, p.y - r].into(), uv: WHITE_UV, color };
        let v2 = Vertex { pos: [p.x + r, p.y + r].into(), uv: WHITE_UV, color };
        let v3 = Vertex { pos: [p.x - r, p.y + r].into(), uv: WHITE_UV, color };
        mesh.vertices.extend_from_slice(&[v0, v1, v2, v3]);
        mesh.indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
    }

    /// Build subduction point meshes (one mesh per color), subsampled to a cap.
    pub fn rebuild_subduction_meshes(
        &mut self,
        rect: Rect,
        latlon: &[[f32; 2]],
        masks: &engine::subduction::SubductionMasks,
    ) {
        let cap = self.effective_subd_cap().max(1) as usize;
        let mut idx_trench: Vec<usize> = Vec::new();
        let mut idx_arc: Vec<usize> = Vec::new();
        let mut idx_back: Vec<usize> = Vec::new();
        for i in 0..latlon.len() {
            if i >= masks.trench.len() {
                break;
            }
            if masks.trench[i] {
                idx_trench.push(i);
            } else if masks.arc[i] {
                idx_arc.push(i);
            } else if masks.backarc[i] {
                idx_back.push(i);
            }
        }
        let total = idx_trench.len() + idx_arc.len() + idx_back.len();
        let total = total.max(1);
        let budget_t =
            ((cap as f32) * (idx_trench.len() as f32) / (total as f32)).round().max(1.0) as usize;
        let budget_a =
            ((cap as f32) * (idx_arc.len() as f32) / (total as f32)).round().max(1.0) as usize;
        let budget_b =
            ((cap as f32) * (idx_back.len() as f32) / (total as f32)).round().max(1.0) as usize;
        let stride_t = (idx_trench.len() / budget_t.max(1)).max(1);
        let stride_a = (idx_arc.len() / budget_a.max(1)).max(1);
        let stride_b = (idx_back.len() / budget_b.max(1)).max(1);

        let mut mesh_t = Mesh::default();
        let mut mesh_a = Mesh::default();
        let mut mesh_b = Mesh::default();
        let col_t = Color32::from_rgb(255, 0, 255); // magenta
        let col_a = Color32::from_rgb(255, 165, 0); // orange
        let col_b = Color32::from_rgb(0, 200, 0); // green
        let r = 1.5;

        for &i in idx_trench.iter().step_by(stride_t) {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_t, p, r, col_t);
        }
        for &i in idx_arc.iter().step_by(stride_a) {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_a, p, r, col_a);
        }
        for &i in idx_back.iter().step_by(stride_b) {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_b, p, r, col_b);
        }

        let mut shapes_t = Vec::new();
        let mut shapes_a = Vec::new();
        let mut shapes_b = Vec::new();
        if !mesh_t.vertices.is_empty() {
            shapes_t.push(Shape::mesh(mesh_t));
        }
        if !mesh_a.vertices.is_empty() {
            shapes_a.push(Shape::mesh(mesh_a));
        }
        if !mesh_b.vertices.is_empty() {
            shapes_b.push(Shape::mesh(mesh_b));
        }
        self.subd_trench = Some(shapes_t);
        self.subd_arc = Some(shapes_a);
        self.subd_backarc = Some(shapes_b);
    }

    /// Build transform point meshes (cyan for pull-apart, brown for restraining).
    pub fn rebuild_transform_meshes(
        &mut self,
        rect: Rect,
        latlon: &[[f32; 2]],
        masks: &engine::transforms::TransformMasks,
    ) {
        let cap = self.trans_max_points.max(1) as usize;
        let mut idx_pull: Vec<usize> = Vec::new();
        let mut idx_rest: Vec<usize> = Vec::new();
        for i in 0..latlon.len() {
            if i >= masks.pull_apart.len() {
                break;
            }
            if masks.pull_apart[i] {
                idx_pull.push(i);
            }
            if masks.restraining[i] {
                idx_rest.push(i);
            }
        }
        let total = (idx_pull.len() + idx_rest.len()).max(1);
        let budget_p =
            ((cap as f32) * (idx_pull.len() as f32) / (total as f32)).round().max(1.0) as usize;
        let budget_r =
            ((cap as f32) * (idx_rest.len() as f32) / (total as f32)).round().max(1.0) as usize;
        let stride_p = (idx_pull.len() / budget_p.max(1)).max(1);
        let stride_r = (idx_rest.len() / budget_r.max(1)).max(1);

        let mut mesh_p = Mesh::default();
        let mut mesh_r = Mesh::default();
        let col_p = Color32::from_rgb(0, 255, 255); // cyan
        let col_r = Color32::from_rgb(150, 75, 0); // brown
        let r = 1.5;
        for &i in idx_pull.iter().step_by(stride_p) {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_p, p, r, col_p);
        }
        for &i in idx_rest.iter().step_by(stride_r) {
            let p = project_equirect(latlon[i][0], latlon[i][1], rect);
            Self::mesh_add_dot(&mut mesh_r, p, r, col_r);
        }
        let mut shapes_p = Vec::new();
        let mut shapes_r = Vec::new();
        if !mesh_p.vertices.is_empty() {
            shapes_p.push(Shape::mesh(mesh_p));
        }
        if !mesh_r.vertices.is_empty() {
            shapes_r.push(Shape::mesh(mesh_r));
        }
        self.trans_pull = Some(shapes_p);
        self.trans_rest = Some(shapes_r);
    }
}
