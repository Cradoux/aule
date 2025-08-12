//! Aulë viewer binary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

mod overlay;

use egui_wgpu::Renderer as EguiRenderer;
use egui_wgpu::ScreenDescriptor;
use egui_winit::State as EguiWinitState;
use winit::{
    dpi::PhysicalSize,
    event::{Event, WindowEvent},
    event_loop::EventLoop,
    window::{Window, WindowBuilder},
};
fn log_grid_info() {
    let f: u32 = 64;
    // Build or load cache (path-agnostic in engine; here we just build).
    let g = engine::grid::Grid::new(f);
    let tiling = engine::grid::tile::Tiling::new(&g, 2, 8192);
    let cells = g.cells as u32; // expected 10*F^2+2
    let pent = 12u32;
    let hex = cells.saturating_sub(pent);
    let mean: f64 = g.area.iter().map(|a| *a as f64).sum::<f64>() / g.area.len() as f64;
    let mut areas: Vec<f32> = g.area.clone();
    areas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if areas.is_empty() { 0.0 } else { areas[areas.len() / 2] } as f64;
    let n1_sample =
        if !g.n1.is_empty() { format!("{:?}", g.n1[0].as_slice()) } else { "[]".to_string() };
    let tiles = tiling.tiles.len();
    let mut inter_sizes: Vec<usize> = tiling.tiles.iter().map(|t| t.interior.len()).collect();
    inter_sizes.sort_unstable();
    let min_i = inter_sizes.first().copied().unwrap_or(0);
    let max_i = inter_sizes.last().copied().unwrap_or(0);
    let mean_i = if inter_sizes.is_empty() {
        0.0
    } else {
        inter_sizes.iter().sum::<usize>() as f64 / inter_sizes.len() as f64
    };
    let mut halo_sizes: Vec<usize> = tiling.tiles.iter().map(|t| t.halo.len()).collect();
    halo_sizes.sort_unstable();
    let min_h = halo_sizes.first().copied().unwrap_or(0);
    let max_h = halo_sizes.last().copied().unwrap_or(0);
    let mean_h = if halo_sizes.is_empty() {
        0.0
    } else {
        halo_sizes.iter().sum::<usize>() as f64 / halo_sizes.len() as f64
    };
    let neigh0 = if !tiling.tiles.is_empty() {
        format!("{:?}", tiling.tiles[0].neighbors.as_slice())
    } else {
        "[]".to_string()
    };
    println!(
        "[grid] F={} cells={} pentagons={} hexagons={} mean_area={:.6} median_area={:.6} n1[0]={} | tiles={} interior[min/mean/max]=[{}/{:.1}/{}] halo[min/mean/max]=[{}/{:.1}/{}] neighbors(tile0)={}",
        f, cells, pent, hex, mean, median, n1_sample, tiles, min_i, mean_i, max_i, min_h, mean_h, max_h, neigh0
    );
}

struct GpuState<'w> {
    _instance: wgpu::Instance,
    surface: wgpu::Surface<'w>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
}
impl<'w> GpuState<'w> {
    #[allow(dead_code)]
    async fn new(window: &'w Window) -> Self {
        let size = window.inner_size();
        let instance = wgpu::Instance::default();
        let surface = match instance.create_surface(window) {
            Ok(s) => s,
            Err(e) => panic!("create surface: {e}"),
        };

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .unwrap_or_else(|| panic!("no suitable GPU adapters"));

        let required_limits = wgpu::Limits::default();

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("device"),
                    required_features: wgpu::Features::empty(),
                    required_limits,
                },
                None,
            )
            .await
            .unwrap_or_else(|e| panic!("request device: {e}"));

        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(surface_caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        Self { _instance: instance, surface, device, queue, config }
    }

    fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
        }
    }

    #[allow(dead_code)]
    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let frame = self.surface.get_current_texture()?;
        let view = frame.texture.create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("encoder") });

        {
            let _rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("clear pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.02,
                            g: 0.02,
                            b: 0.04,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
        }

        self.queue.submit(std::iter::once(encoder.finish()));
        frame.present();

        Ok(())
    }
}
fn main() {
    let event_loop = EventLoop::new().unwrap_or_else(|e| panic!("event loop: {e}"));
    let title = format!("Aulë Viewer v{}", engine::version());
    let window_init = WindowBuilder::new()
        .with_title(title)
        .build(&event_loop)
        .unwrap_or_else(|e| panic!("create window: {e}"));

    // Leak the window to obtain a 'static reference for the surface lifetime without unsafe.
    let window: &'static Window = Box::leak(Box::new(window_init));
    let mut gpu = pollster::block_on(GpuState::new(window));
    let egui_ctx = egui::Context::default();
    let mut egui_state =
        EguiWinitState::new(egui_ctx.clone(), egui::ViewportId::ROOT, &event_loop, None, None);
    let surface_format = gpu.config.format;
    let mut egui_renderer = EguiRenderer::new(&gpu.device, surface_format, None, 1);
    let mut ov = overlay::OverlayState::default();
    // T-020: Construct device field buffers sized to the grid (then drop)
    {
        let f: u32 = 64;
        let g_tmp = engine::grid::Grid::new(f);
        let _device_fields = engine::fields::DeviceFields::new(&gpu.device, g_tmp.cells);
    }
    // Build data for overlays once
    let f: u32 = 64;
    let g_view = engine::grid::Grid::new(f);
    let plates = engine::plates::Plates::new(&g_view, 8, 12345);
    let mut mags: Vec<f64> =
        plates.vel_en.iter().map(|v| ((v[0] as f64).hypot(v[1] as f64))).collect();
    mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = mags.len();
    let min_v = *mags.first().unwrap_or(&0.0);
    let max_v = *mags.last().unwrap_or(&0.0);
    let mean_v = if n == 0 { 0.0 } else { mags.iter().sum::<f64>() / n as f64 };
    println!(
        "[plates] N={} |V| min/mean/max = {:.3} / {:.3} / {:.3} m/yr",
        8, min_v, mean_v, max_v
    );
    let bounds =
        engine::boundaries::Boundaries::classify(&g_view, &plates.plate_id, &plates.vel_en, 0.005);
    println!(
        "[boundaries] div={} conv={} trans={} (τ=0.5 cm/yr)",
        bounds.stats.divergent, bounds.stats.convergent, bounds.stats.transform
    );

    // Ridge CPU pass: initialize a host age buffer and apply births + fringe
    let mut age_ocean = vec![10.0f32; g_view.cells];
    let ridge_stats = engine::ridge::apply_ridge(
        &g_view,
        &bounds,
        &mut age_ocean,
        engine::ridge::RidgeParams { fringe_age_myr: 0.2 },
    );
    println!(
        "[ridge] births={} fringe={} (fringe_age={} Myr)",
        ridge_stats.births, ridge_stats.fringe, 0.2
    );

    // Compute steady-state age/bathymetry once (CPU)
    let age_params =
        engine::age::AgeParams { v_floor_m_per_yr: (ov.v_floor_cm_per_yr as f64) * 0.01 };
    let age_out = engine::age::compute_age_and_bathymetry(
        &g_view,
        &bounds,
        &plates.plate_id,
        &plates.vel_en,
        age_params,
    );
    ov.age_minmax = age_out.min_max_age;
    ov.depth_minmax = age_out.min_max_depth;

    log_grid_info();

    let mut last_frame = std::time::Instant::now();
    let mut fps: f32 = 0.0;
    let nplates: usize = plates.pole_axis.len();

    event_loop
    .run(move |event, elwt| {
        match event {
            Event::AboutToWait => {
                window.request_redraw();
            }
            Event::WindowEvent { event, window_id } if window_id == window.id() => {
                // forward events to egui (note: window, not context)
                let _ = egui_state.on_window_event(window, &event);
                match event {
                    WindowEvent::CloseRequested => elwt.exit(),
                    WindowEvent::Resized(size) => {
                        gpu.resize(size);
                    }
                    WindowEvent::RedrawRequested => {
                        let raw_input = egui_state.take_egui_input(window);
                        let full_output = egui_ctx.run(raw_input, |ctx| {
                            if ctx.input(|i| i.key_pressed(egui::Key::Num1)) { ov.show_plates = !ov.show_plates; ov.plates_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num2)) { ov.show_vel = !ov.show_vel; ov.vel_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num3)) { ov.show_bounds = !ov.show_bounds; ov.bounds_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num4)) { ov.show_age = !ov.show_age; ov.age_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num5)) { ov.show_bathy = !ov.show_bathy; ov.bathy_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::H)) { ov.show_hud = !ov.show_hud; }

                            egui::TopBottomPanel::top("hud").show_animated(ctx, ov.show_hud, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("1: Plates  2: Velocities  3: Boundaries  4: Age  5: Bathy  H: HUD");
                                    ui.separator();
                                    ui.label(format!(
                                        "plates={}  |V| min/mean/max = {:.2}/{:.2}/{:.2} cm/yr",
                                        nplates,
                                        min_v * 100.0,
                                        mean_v * 100.0,
                                        max_v * 100.0
                                    ));
                                    ui.separator();
                                    ui.label(format!(
                                        "boundaries: div={} conv={} trans={}",
                                        bounds.stats.divergent, bounds.stats.convergent, bounds.stats.transform
                                    ));
                                    ui.separator();
                                    ui.label(format!(
                                        "max arrows={} max strokes={}  scale={:.2} px/cm/yr  v_floor={:.2} cm/yr  FPS: {:.0}",
                                        ov.max_arrows_slider, ov.max_bounds_slider, ov.vel_scale_px_per_cm_yr, ov.v_floor_cm_per_yr, fps
                                    ));
                                });
                                ui.separator();
                                ui.horizontal(|ui| {
                                    ui.add(
                                        egui::Slider::new(&mut ov.vel_scale_px_per_cm_yr, 0.1..=2.0)
                                            .text("Vel scale (px per cm/yr)")
                                    );
                                    let arrows = egui::Slider::new(&mut ov.max_arrows_slider, 500..=20_000)
                                        .text("Max arrows").step_by(500.0);
                                    ui.add(arrows);
                                    let bounds_cap = egui::Slider::new(&mut ov.max_bounds_slider, 500..=20_000)
                                        .text("Max boundaries").step_by(500.0);
                                    ui.add(bounds_cap);
                                });
                                ui.horizontal(|ui| {
                                    ui.checkbox(&mut ov.adaptive_cap, "Adaptive cap (16.6 ms target)");
                                    ui.label(format!(
                                        "live arrows={} live boundaries={}",
                                        ov.live_arrows_cap, ov.live_bounds_cap
                                    ));
                                });
                                ui.separator();
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(format!(
                                        "age min/mean/max = {:.2}/{:.2}/{:.2} Myr",
                                        age_out.min_max_age.0,
                                        (age_out.age_myr.iter().sum::<f32>() / age_out.age_myr.len().max(1) as f32),
                                        age_out.min_max_age.1
                                    ));
                                    ui.separator();
                                    ui.label(format!(
                                        "depth min/mean/max = {:.0}/{:.0}/{:.0} m",
                                        age_out.min_max_depth.0,
                                        (age_out.depth_m.iter().sum::<f32>() / age_out.depth_m.len().max(1) as f32),
                                        age_out.min_max_depth.1
                                    ));
                                });
                            });

                            egui::CentralPanel::default().show(ctx, |ui| {
                                let rect = ui.available_rect_before_wrap();
                                let painter = ui.painter_at(rect);
                                // Ensure caches are valid for current params
                                let eff_ar = ov.effective_arrows_cap();
                                let eff_bd = ov.effective_bounds_cap();
                                ov.ensure_params_and_invalidate_if_needed(rect, eff_ar, eff_bd);
                                if ov.show_plates {
                                    for s in ov.shapes_for_plates(rect, &g_view.latlon, &plates.plate_id) { painter.add(s.clone()); }
                                }
                                if ov.show_vel {
                                    for s in ov.shapes_for_velocities(rect, &g_view.latlon, &plates.vel_en) { painter.add(s.clone()); }
                                }
                                if ov.show_bounds {
                                    for s in ov.shapes_for_boundaries(rect, &g_view.latlon, &bounds.edges) { painter.add(s.clone()); }
                                }
                                if ov.show_age {
                                    if ov.age_cache.is_none() { ov.rebuild_age_shapes(rect, &g_view.latlon, &age_out.age_myr); }
                                    for s in ov.age_shapes() { painter.add(s.clone()); }
                                }
                                if ov.show_bathy {
                                    if ov.bathy_cache.is_none() { ov.rebuild_bathy_shapes(rect, &g_view.latlon, &age_out.depth_m); }
                                    for s in ov.bathy_shapes() { painter.add(s.clone()); }
                                }
                            });
                        });

                        for (id, image_delta) in &full_output.textures_delta.set {
                            egui_renderer.update_texture(&gpu.device, &gpu.queue, *id, image_delta);
                        }
                        for id in &full_output.textures_delta.free {
                            egui_renderer.free_texture(id);
                        }
                        let ppp = window.scale_factor() as f32;
                        let paint_jobs = egui_ctx.tessellate(full_output.shapes, ppp);

                        let frame = match gpu.surface.get_current_texture() {
                            Ok(f) => f,
                            Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => { gpu.resize(window.inner_size()); return; }
                            Err(wgpu::SurfaceError::OutOfMemory) => { elwt.exit(); return; }
                            Err(wgpu::SurfaceError::Timeout) => { return; }
                        };
                        let view = frame.texture.create_view(&wgpu::TextureViewDescriptor::default());
                        let mut encoder = gpu
                            .device
                            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("encoder") });

                        let screen_desc = ScreenDescriptor {
                            size_in_pixels: [gpu.config.width, gpu.config.height],
                            pixels_per_point: ppp,
                        };
                        egui_renderer.update_buffers(
                            &gpu.device,
                            &gpu.queue,
                            &mut encoder,
                            &paint_jobs,
                            &screen_desc,
                        );

                        {
                            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                                label: Some("egui pass"),
                                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                                    view: &view,
                                    resolve_target: None,
                                    ops: wgpu::Operations {
                                        load: wgpu::LoadOp::Load,
                                        store: wgpu::StoreOp::Store,
                                    },
                                })],
                                depth_stencil_attachment: None,
                                occlusion_query_set: None,
                                timestamp_writes: None,
                            });
                            egui_renderer.render(&mut rpass, &paint_jobs, &screen_desc);
                        }
                        gpu.queue.submit(std::iter::once(encoder.finish()));
                        frame.present();

                        egui_state.handle_platform_output(window, full_output.platform_output);
                        let now = std::time::Instant::now();
                        let dt = now.duration_since(last_frame).as_secs_f32();
                        last_frame = now;
                        if dt > 0.0 { fps = 0.9 * fps + 0.1 * (1.0 / dt); }
                        // Update adaptive caps using frame time in ms
                        ov.update_adaptive_caps(dt * 1000.0);
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    })
    .unwrap_or_else(|e| panic!("run app: {e}"));
}
