//! Aulë viewer binary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

mod overlay;
mod plot;

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
    // Build world state and initial overlay data
    let f: u32 = 64;
    let mut world = engine::world::World::new(f, 8, 12345);
    let mut mags: Vec<f64> =
        world.v_en.iter().map(|v| ((v[0] as f64).hypot(v[1] as f64))).collect();
    mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = mags.len();
    let min_v = *mags.first().unwrap_or(&0.0);
    let max_v = *mags.last().unwrap_or(&0.0);
    let mean_v = if n == 0 { 0.0 } else { mags.iter().sum::<f64>() / n as f64 };
    println!(
        "[plates] N={} |V| min/mean/max = {:.3} / {:.3} / {:.3} m/yr",
        8, min_v, mean_v, max_v
    );
    println!(
        "[boundaries] div={} conv={} trans={} (τ=0.5 cm/yr)",
        world.boundaries.stats.divergent,
        world.boundaries.stats.convergent,
        world.boundaries.stats.transform
    );
    // Initialize overlay ranges
    ov.age_minmax = (0.0, 100.0);
    ov.depth_minmax = (2600.0, 6000.0);

    log_grid_info();

    let mut last_frame = std::time::Instant::now();
    let mut fps: f32 = 0.0;
    let nplates: usize = world.plates.pole_axis.len();

    // Playback state
    let mut playing: bool = false;
    let mut dt_myr: f32 = 1.0;
    let mut steps_per_sec: u32 = 5;
    let mut step_accum_s: f32 = 0.0;

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
                            if ctx.input(|i| i.key_pressed(egui::Key::Num6)) { ov.show_age_depth = !ov.show_age_depth; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num7)) { ov.show_subduction = !ov.show_subduction; ov.subd_trench=None; ov.subd_arc=None; ov.subd_backarc=None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::Num0)) { ov.show_transforms = !ov.show_transforms; ov.trans_pull=None; ov.trans_rest=None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::H)) { ov.show_hud = !ov.show_hud; }

                            egui::TopBottomPanel::top("hud").show_animated(ctx, ov.show_hud, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("1: Plates  2: Velocities  3: Boundaries  4: Age  5: Bathy  6: Age–Depth  7: Subduction  H: HUD");
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
                                        "boundaries: div={} conv={} trans={} | transforms: pull_apart={} restraining={}",
                                        world.boundaries.stats.divergent, world.boundaries.stats.convergent, world.boundaries.stats.transform,
                                        ov.trans_pull_count, ov.trans_rest_count
                                    ));
                                    ui.separator();
                                    if playing { if ui.button("⏸").clicked() { playing = false; } } else if ui.button("▶").clicked() { playing = true; }
                                    if ui.button("⏭").clicked() {
                                        let p = engine::stepper::StepParams { dt_myr, tau_open_m_per_yr: 0.005 };
                                        let _ = engine::stepper::step(&mut world, &p);
                                        ov.age_cache=None; ov.bathy_cache=None; ov.bounds_cache=None; ov.subd_trench=None; ov.subd_arc=None; ov.subd_backarc=None;
                                    }
                                    ui.separator();
                                    ui.label(format!(
                                        "max arrows={} max strokes={} max subd={}  scale={:.2} px/cm/yr  v_floor={:.2} cm/yr  FPS: {:.0}",
                                        ov.max_arrows_slider, ov.max_bounds_slider, ov.max_subd_slider, ov.vel_scale_px_per_cm_yr, ov.v_floor_cm_per_yr, fps
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
                                    let subd_cap = egui::Slider::new(&mut ov.max_subd_slider, 500..=20_000)
                                        .text("Max subduction points").step_by(500.0);
                                    ui.add(subd_cap);
                                    ui.separator();
                                    ui.add(egui::Slider::new(&mut dt_myr, 0.1..=5.0).text("dt (Myr)"));
                                    ui.add(egui::Slider::new(&mut steps_per_sec, 1..=30).text("speed (steps/sec)"));
                                    ui.label(format!("t={:.1} Myr  step={}", world.clock.t_myr, world.clock.step_idx));
                                });
                                ui.separator();
                                ui.collapsing("Transforms (0)", |ui| {
                                    ui.label("Active if |tangential| ≥ min_tangential & |normal| ≤ τ_open");
                                    ui.label("Cyan = pull-apart (deeper), Brown = restraining (shallower)");
                                    ui.label("Width = half-width from the fault trace (km)");
                                    let mut changed = false;
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_min_tangential_m_per_yr, 0.0001..=0.03).logarithmic(true).text("min tangential (m/yr)")).changed();
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_tau_open_m_per_yr, 0.002..=0.02).logarithmic(true).text("τ_open (m/yr)")).changed();
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_basin_half_width_km, 10.0..=60.0).text("Half-width (km)")).changed();
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_basin_deepen_m, 100.0..=1000.0).text("Basin deepen (m)")).changed();
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_ridge_like_uplift_m, -800.0..=-50.0).text("Restraining uplift (m)")).changed();
                                    changed |= ui.add(egui::Slider::new(&mut ov.trans_max_points, 500..=20_000).text("Max points")).changed();
                                    if changed {
                                        // Invalidate caches so we recompute once below
                                        ov.trans_pull=None; ov.trans_rest=None; ov.bathy_cache=None;
                                    }
                                });
                                ui.horizontal(|ui| {
                                    ui.checkbox(&mut ov.adaptive_cap, "Adaptive cap (16.6 ms target)");
                                    ui.label(format!(
                                        "live arrows={} live boundaries={} live subd={}",
                                        ov.live_arrows_cap, ov.live_bounds_cap, ov.live_subd_cap
                                    ));
                                });
                                ui.separator();
                                ui.horizontal_wrapped(|ui| {
                                    let (mut amin, mut amax) = (f32::INFINITY, f32::NEG_INFINITY);
                                    let (mut dmin, mut dmax) = (f32::INFINITY, f32::NEG_INFINITY);
                                    for &a in &world.age_myr {
                                        if a.is_finite() {
                                            if a < amin { amin = a; }
                                            if a > amax { amax = a; }
                                        }
                                    }
                                    for &d in &world.depth_m {
                                        if d.is_finite() {
                                            if d < dmin { dmin = d; }
                                            if d > dmax { dmax = d; }
                                        }
                                    }
                                    let amin = if amin.is_finite() { amin } else { 0.0 };
                                    let amax = if amax.is_finite() { amax } else { 0.0 };
                                    let dmin = if dmin.is_finite() { dmin } else { 0.0 };
                                    let dmax = if dmax.is_finite() { dmax } else { 0.0 };
                                    let amean = world.age_myr.iter().copied().sum::<f32>() / world.age_myr.len().max(1) as f32;
                                    let dmean = world.depth_m.iter().copied().sum::<f32>() / world.depth_m.len().max(1) as f32;
                                    ui.label(format!("age min/mean/max = {:.2}/{:.2}/{:.2} Myr", amin, amean, amax));
                                    ui.separator();
                                    ui.label(format!("depth min/mean/max = {:.0}/{:.0}/{:.0} m", dmin, dmean, dmax));
                                });
                            });

                            // Optional bottom panel for age-depth plot
                            if ov.show_age_depth {
                                egui::TopBottomPanel::bottom("age_depth_panel").show(ctx, |ui| {
                                    ui.label("Age–Depth Validation");
                                    ui.horizontal(|ui| {
                                        ui.add(
                                            egui::Slider::new(&mut ov.plot_sample_cap, 1000..=50_000)
                                                .text("Sample cap")
                                                .step_by(1000.0),
                                        );
                                        ui.add(
                                            egui::Slider::new(&mut ov.plot_bin_width_myr, 1.0..=20.0)
                                                .text("Bin width (Myr)"),
                                        );
                                    });
                                    let d0 = 2600.0_f64;
                                    let a = 350.0_f64;
                                    let b = 0.0_f64;
                                    let pdata = plot::build_age_depth_plot(
                                        &world.age_myr,
                                        &world.depth_m,
                                        plot::AgeDepthPlotParams { sample_cap: ov.plot_sample_cap as usize, bin_width_myr: ov.plot_bin_width_myr },
                                        &|age| engine::age::depth_from_age(age, d0, a, b),
                                    );
                                    ui.label("Sample cap = #points plotted (subsampled).\nBin width = age range per bin for the yellow binned mean curve; smaller = more detailed/noisier.");
                                    ui.horizontal(|ui| {
                                        ui.label(format!("RMS: {:.2} m  N={}  age[min/max]={:.2}/{:.2} Myr  depth[min/max]={:.0}/{:.0} m",
                                            pdata.stats.rms_m, pdata.stats.n_samples,
                                            pdata.stats.age_minmax.0, pdata.stats.age_minmax.1,
                                            pdata.stats.depth_minmax.0, pdata.stats.depth_minmax.1));
                                    });
                                    egui_plot::Plot::new("age_depth_plot")
                                        .legend(egui_plot::Legend::default())
                                        .x_axis_label("Age (Myr)")
                                        .y_axis_label("Depth (m)")
                                        .show(ui, |plot_ui| {
                                            plot_ui.points(pdata.scatter);
                                            plot_ui.line(pdata.binned);
                                            plot_ui.line(pdata.reference);
                                        });
                                });
                            }

                            egui::CentralPanel::default().show(ctx, |ui| {
                                let rect = ui.available_rect_before_wrap();
                                let painter = ui.painter_at(rect);
                                // Ensure caches are valid for current params
                                let eff_ar = ov.effective_arrows_cap();
                                let eff_bd = ov.effective_bounds_cap();
                                let eff_sd = ov.effective_subd_cap();
                                ov.ensure_params_and_invalidate_if_needed(rect, eff_ar, eff_bd, eff_sd);
                                if ov.show_plates {
                                    for s in ov.shapes_for_plates(rect, &world.grid.latlon, &world.plates.plate_id) { painter.add(s.clone()); }
                                }
                                if ov.show_vel {
                                    for s in ov.shapes_for_velocities(rect, &world.grid.latlon, &world.v_en) { painter.add(s.clone()); }
                                }
                                if ov.show_bounds {
                                    for s in ov.shapes_for_boundaries(rect, &world.grid.latlon, &world.boundaries.edges) { painter.add(s.clone()); }
                                }
                                if ov.show_age {
                                    if ov.age_cache.is_none() { ov.rebuild_age_shapes(rect, &world.grid.latlon, &world.age_myr); }
                                    for s in ov.age_shapes() { painter.add(s.clone()); }
                                }
                                if ov.show_bathy {
                                    if ov.bathy_cache.is_none() { ov.rebuild_bathy_shapes(rect, &world.grid.latlon, &world.depth_m); }
                                    for s in ov.bathy_shapes() { painter.add(s.clone()); }
                                }
                                if ov.show_subduction {
                                    if ov.subd_trench.is_none() && ov.subd_arc.is_none() && ov.subd_backarc.is_none() {
                                        let mut tmp_depth = vec![0.0f32; world.grid.cells];
                                        let sub_res = engine::subduction::apply_subduction(
                                            &world.grid,
                                            &world.boundaries,
                                            &world.plates.plate_id,
                                            &world.age_myr,
                                            &world.v_en,
                                            &mut tmp_depth,
                                            engine::subduction::SubductionParams {
                                                tau_conv_m_per_yr: 0.005,
                                                trench_half_width_km: 50.0,
                                                arc_offset_km: 150.0,
                                                arc_half_width_km: 30.0,
                                                backarc_width_km: 150.0,
                                                trench_deepen_m: 3000.0,
                                                arc_uplift_m: -500.0,
                                                backarc_uplift_m: -200.0,
                                                rollback_offset_m: 0.0,
                                                rollback_rate_km_per_myr: 0.0,
                                                backarc_extension_mode: false,
                                                backarc_extension_deepen_m: 600.0,
                                            },
                                        );
                                        ov.rebuild_subduction_meshes(rect, &world.grid.latlon, &sub_res.masks);
                                    }
                                    if let Some(v) = &ov.subd_trench { for s in v { painter.add(s.clone()); } }
                                    if let Some(v) = &ov.subd_arc { for s in v { painter.add(s.clone()); } }
                                    if let Some(v) = &ov.subd_backarc { for s in v { painter.add(s.clone()); } }
                                }
                                if ov.show_transforms {
                                    if ov.trans_pull.is_none() && ov.trans_rest.is_none() {
                                        // Recompute depth: baseline -> subduction -> transforms, avoiding overlaps
                                        // Baseline from age
                                        let mut depth_base = vec![0.0f32; world.grid.cells];
                                        for (i, db) in depth_base.iter_mut().enumerate().take(world.grid.cells) {
                                            let mut d = engine::age::depth_from_age(world.age_myr[i] as f64, 2600.0, 350.0, 0.0) as f32;
                                            if !d.is_finite() { d = 6000.0; }
                                            *db = d.clamp(0.0, 6000.0);
                                        }
                                        // Subduction on baseline
                                        let mut depth_subd = depth_base.clone();
                                        let sub_res = engine::subduction::apply_subduction(
                                            &world.grid, &world.boundaries, &world.plates.plate_id,
                                            &world.age_myr, &world.v_en, &mut depth_subd,
                                            engine::subduction::SubductionParams {
                                                tau_conv_m_per_yr: 0.005,
                                                trench_half_width_km: 50.0,
                                                arc_offset_km: 150.0,
                                                arc_half_width_km: 30.0,
                                                backarc_width_km: 150.0,
                                                trench_deepen_m: 3000.0,
                                                arc_uplift_m: -500.0,
                                                backarc_uplift_m: -200.0,
                                                rollback_offset_m: 0.0,
                                                rollback_rate_km_per_myr: 0.0,
                                                backarc_extension_mode: false,
                                                backarc_extension_deepen_m: 600.0,
                                            },
                                        );
                                        // Transforms: compute delta from baseline
                                        let mut depth_trans = depth_base.clone();
                                        let (trans_masks, trans_stats) = engine::transforms::apply_transforms(
                                            &world.grid, &world.boundaries, &world.plates.plate_id, &world.v_en, &mut depth_trans,
                                            engine::transforms::TransformParams {
                                                tau_open_m_per_yr: ov.trans_tau_open_m_per_yr as f64,
                                                min_tangential_m_per_yr: ov.trans_min_tangential_m_per_yr as f64,
                                                basin_half_width_km: ov.trans_basin_half_width_km as f64,
                                                ridge_like_uplift_m: ov.trans_ridge_like_uplift_m,
                                                basin_deepen_m: ov.trans_basin_deepen_m,
                                            }
                                        );
                                        println!("[transforms] bands: pull_apart={} restraining={}", trans_stats.pull_apart_cells, trans_stats.restraining_cells);
                                        ov.rebuild_transform_meshes(rect, &world.grid.latlon, &trans_masks);
                                        ov.trans_pull_count = trans_stats.pull_apart_cells;
                                        ov.trans_rest_count = trans_stats.restraining_cells;
                                        let mut final_depth = depth_subd;
                                        for i in 0..world.grid.cells {
                                            let delta_t = depth_trans[i] - depth_base[i];
                                            // avoid overlap with subduction masks
                                            if !sub_res.masks.trench[i] && !sub_res.masks.arc[i] && !sub_res.masks.backarc[i] {
                                                final_depth[i] += delta_t;
                                            }
                                        }
                                        world.depth_m = final_depth;
                                        ov.bathy_cache = None;
                                    }
                                    if let Some(v) = &ov.trans_pull { for s in v { painter.add(s.clone()); } }
                                    if let Some(v) = &ov.trans_rest { for s in v { painter.add(s.clone()); } }
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
                        // Playback ticking
                        if playing {
                            step_accum_s += dt;
                            let step_interval = 1.0f32 / (steps_per_sec.max(1) as f32);
                            let mut steps_this_frame = 0u32;
                            let max_steps_frame = 4u32; // clamp to keep FPS ~60
                            while step_accum_s >= step_interval && steps_this_frame < max_steps_frame {
                                let p = engine::stepper::StepParams { dt_myr, tau_open_m_per_yr: 0.005 };
                                let _ = engine::stepper::step(&mut world, &p);
                                step_accum_s -= step_interval;
                                steps_this_frame += 1;
                            }
                            if steps_this_frame > 0 {
                                ov.age_cache=None; ov.bathy_cache=None; ov.bounds_cache=None; ov.subd_trench=None; ov.subd_arc=None; ov.subd_backarc=None;
                            }
                        }
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    })
    .unwrap_or_else(|e| panic!("run app: {e}"));
}
