pub struct OrbitCamera {
    pub yaw: f32,
    pub pitch: f32,
    pub distance: f32,
    pub fov_y: f32,
    pub aspect: f32,
    pub z_near: f32,
    pub z_far: f32,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            yaw: 0.6,
            pitch: 0.4,
            distance: 3.0,
            fov_y: 60f32.to_radians(),
            aspect: 1.6,
            z_near: 0.01,
            z_far: 100.0,
        }
    }
}

impl OrbitCamera {
    pub fn update_from_input(&mut self, ctx: &egui::Context, ui_hijacked: bool) {
        if ui_hijacked {
            return;
        }
        ctx.input(|i| {
            let dragging_rmb = i.pointer.button_down(egui::PointerButton::Secondary)
                || (i.modifiers.alt && i.pointer.button_down(egui::PointerButton::Primary));
            if dragging_rmb {
                let d = i.pointer.delta();
                let k = 0.005f32;
                self.yaw -= d.x * k;
                self.pitch -= d.y * k;
                let lim = core::f32::consts::FRAC_PI_2 - 0.017;
                if self.pitch > lim {
                    self.pitch = lim;
                }
                if self.pitch < -lim {
                    self.pitch = -lim;
                }
            }
            let scroll = i.smooth_scroll_delta.y + i.raw_scroll_delta.y;
            if scroll.abs() > 0.0 {
                let factor = (-scroll * 0.0015).exp();
                self.distance = (self.distance * factor).clamp(0.8, 5.0);
            }
        });
    }
    pub fn view_proj(&self) -> [[f32; 4]; 4] {
        let eye = glam::Vec3::new(
            self.distance * self.yaw.cos() * self.pitch.cos(),
            self.distance * self.pitch.sin(),
            self.distance * self.yaw.sin() * self.pitch.cos(),
        );
        let center = glam::Vec3::ZERO;
        let up = glam::Vec3::Y;
        let view = glam::Mat4::look_at_rh(eye, center, up);
        let proj =
            glam::Mat4::perspective_rh(self.fov_y, self.aspect.max(1e-3), self.z_near, self.z_far);
        (proj * view).to_cols_array_2d()
    }

    pub fn eye(&self) -> glam::Vec3 {
        glam::Vec3::new(
            self.distance * self.yaw.cos() * self.pitch.cos(),
            self.distance * self.pitch.sin(),
            self.distance * self.yaw.sin() * self.pitch.cos(),
        )
    }
}
