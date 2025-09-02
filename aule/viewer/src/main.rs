//! Aulë viewer binary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

mod colormap;
mod globe;
mod overlay;
mod plot;
mod plot_age_depth;
mod plot_flexure;
mod raster;
mod raster_gpu;
use viewer::cpu_picker::{GeoPicker, Vec3};
use viewer::pixel_map::{pixel_to_lon_lat, sph_to_unit as sph_to_unit_px};

use egui_wgpu::Renderer as EguiRenderer;
use egui_wgpu::ScreenDescriptor;
use egui_winit::State as EguiWinitState;
use winit::{
    dpi::PhysicalSize,
    event::{Event, WindowEvent},
    event_loop::EventLoop,
    window::{Window, WindowBuilder},
};

// T-505 drawer render — Simple/Advanced
#[allow(dead_code)]
fn run_to_t_realtime(
    ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
    sp: &engine::world::StepParams,
    t_end_myr: f64,
    max_steps_per_yield: u32,
) {
    let max_chunk = max_steps_per_yield.max(1);
    while world.clock.t_myr < t_end_myr {
        for _ in 0..max_chunk {
            if world.clock.t_myr >= t_end_myr {
                break;
            }
            // Centralized pipeline step: solves eta to target land fraction
            let cfg = engine::pipeline::PipelineCfg {
                dt_myr: sp.dt_myr as f32,
                steps_per_frame: 1,
                enable_flexure: sp.do_flexure,
                enable_erosion: sp.do_surface,
                target_land_frac: ov.simple_target_land,
                freeze_eta: ov.freeze_eta,
                log_mass_budget: false,
                enable_subduction: sp.do_subduction,
                enable_rigid_motion: sp.do_rigid_motion,
                cadence_trf_every: ov.cadence_trf_every.max(1),
                cadence_sub_every: ov.cadence_sub_every.max(1),
                cadence_flx_every: ov.cadence_flx_every.max(1),
                cadence_sea_every: ov.cadence_sea_every.max(1),
                cadence_surf_every: ov.cadence_sea_every.max(1),
                substeps_transforms: 4,
                substeps_subduction: 4,
                use_gpu_flexure: true,
                gpu_flex_levels: ov.levels.max(1),
                gpu_flex_cycles: ov.flex_cycles.max(1),
                gpu_wj_omega: ov.wj_omega,
                subtract_mean_load: ov.subtract_mean_load,
                surf_k_stream: ov.surf_k_stream,
                surf_m_exp: ov.surf_m_exp,
                surf_n_exp: ov.surf_n_exp,
                surf_k_diff: ov.surf_k_diff,
                surf_k_tr: ov.surf_k_tr,
                surf_p_exp: ov.surf_p_exp,
                surf_q_exp: ov.surf_q_exp,
                surf_rho_sed: ov.surf_rho_sed,
                surf_min_slope: ov.surf_min_slope,
                surf_subcycles: ov.surf_subcycles,
                surf_couple_flexure: ov.surf_couple_flexure,
                sub_tau_conv_m_per_yr: ov.sub_tau_conv_m_per_yr,
                sub_trench_half_width_km: ov.sub_trench_half_width_km,
                sub_arc_offset_km: ov.sub_arc_offset_km,
                sub_arc_half_width_km: ov.sub_arc_half_width_km,
                sub_backarc_width_km: ov.sub_backarc_width_km,
                sub_trench_deepen_m: ov.sub_trench_deepen_m,
                sub_arc_uplift_m: ov.sub_arc_uplift_m,
                sub_backarc_uplift_m: ov.sub_backarc_uplift_m,
                sub_rollback_offset_m: ov.sub_rollback_offset_m,
                sub_rollback_rate_km_per_myr: ov.sub_rollback_rate_km_per_myr,
                sub_backarc_extension_mode: ov.sub_backarc_extension_mode,
                sub_backarc_extension_deepen_m: ov.sub_backarc_extension_deepen_m,
                sub_continent_c_min: ov.sub_continent_c_min,
                cadence_spawn_plate_every: 0,
                cadence_retire_plate_every: 0,
                cadence_force_balance_every: 8,
                fb_gain: 1.0e-12,
                fb_damp_per_myr: 0.2,
                fb_k_conv: 1.0,
                fb_k_div: 0.5,
                fb_k_trans: 0.1,
                fb_max_domega: 5.0e-9,
                fb_max_omega: 2.0e-7,
            };
            let mut elev_tmp: Vec<f32> = vec![0.0; world.depth_m.len()];
            let mut eta = world.sea.eta_m;
            engine::pipeline::step_full(
                world,
                engine::pipeline::SurfaceFields { elev_m: &mut elev_tmp, eta_m: &mut eta },
                cfg,
            );
            world.sea.eta_m = eta;
            ov.world_dirty = true;
            ov.color_dirty = true;
            ov.bathy_cache = None;
            ctx.request_repaint();
        }
        std::thread::yield_now();
    }
}
fn render_simple_panels(
    ui: &mut egui::Ui,
    ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
) {
    egui::CollapsingHeader::new("Simulation Basics").default_open(true).show(ui, |ui| {
        let tiling_ok = cfg!(feature = "tiling") || ov.high_f_available;
        let current_f = world.grid.frequency;
        ui.horizontal(|ui| {
            ui.label("Resolution (F):");
            ui.selectable_value(&mut ov.simple_f, 64, "64");
            ui.selectable_value(&mut ov.simple_f, 128, "128");
            ui.selectable_value(&mut ov.simple_f, 256, "256");
            ui.add_enabled_ui(tiling_ok, |ui| {
                ui.selectable_value(&mut ov.simple_f, 512, "512");
            });
            if !tiling_ok {
                ui.label(egui::RichText::new("512 requires T-455").small());
            }
        });
        // Rebuild grid/world on F change (no auto run)
        if ov.simple_f != current_f {
            let allow = ov.simple_f <= 256 || tiling_ok;
            if allow {
                let num_plates: u32 = world.plates.pole_axis.len() as u32;
                let num_plates = if num_plates == 0 { 8 } else { num_plates };
                *world = engine::world::World::new(ov.simple_f, num_plates, ov.simple_seed);
                ov.world_dirty = true;
                ov.bathy_cache = None;
                ov.raster_dirty = true;
                ctx.request_repaint();
            }
        }
        ui.add(egui::DragValue::new(&mut ov.simple_t_end_myr).speed(10.0).suffix(" Myr"));
        ui.horizontal(|ui| {
            if ui.button("▶ Play").clicked() {
                ov.stepper.playing = true;
                ov.stepper.t_target_myr = ov.simple_t_end_myr as f32;
                ov.run_active = true; // keep resolution policy in sync until fully migrated
                ov.raster_dirty = true;
                ctx.request_repaint();
            }
            if ui.button("⏸ Pause").clicked() {
                ov.stepper.playing = false;
                ov.run_active = false;
                ov.raster_dirty = true;
                ctx.request_repaint();
            }
            ui.label("steps/frame:");
            let mut mspf = ov.stepper.max_steps_per_frame;
            if ui.add(egui::Slider::new(&mut mspf, 1..=8)).changed() {
                ov.stepper.max_steps_per_frame = mspf;
            }
        });

        if ui.button("Generate world").clicked() {
            // Simple-mode generation using area-based land solver (T-603c)
            let plates = match ov.simple_preset {
                1 => 10,
                2 => 6,
                3 => 8, // Supercontinent
                _ => 8,
            };
            let continents_n = match ov.simple_preset {
                1 => 4,
                2 => 2,
                3 => 0, // Supercontinent preset triggers specialized generator
                _ => 3,
            };
            // Reset/normalize world state for current F & plates
            world.plates = engine::plates::Plates::new(&world.grid, plates, ov.simple_seed);
            world.v_en.clone_from(&world.plates.vel_en);
            world.age_myr.fill(0.0);
            world.depth_m.fill(0.0);
            world.c.fill(0.0);
            world.th_c_m.fill(0.0);
            world.sediment_m.fill(0.0);
            world.clock = engine::world::Clock { t_myr: 0.0, step_idx: 0 };
            world.sea_level_ref = None;
            world.boundaries = engine::boundaries::Boundaries::classify(
                &world.grid,
                &world.plates.plate_id,
                &world.v_en,
                0.005,
            );
            // Initialize ridge births then baseline bathymetry from age (ocean depths)
            {
                let mut ages_tmp = world.age_myr.clone();
                let _ridge_stats = engine::ridge::apply_ridge(
                    &world.grid,
                    &world.boundaries,
                    &mut ages_tmp,
                    engine::ridge::RidgeParams { fringe_age_myr: 0.0 },
                );
                world.age_myr = ages_tmp;
            }
            let n_cells = world.grid.cells;
            for i in 0..n_cells {
                let mut d = engine::age::depth_from_age_plate(
                    world.age_myr[i] as f64,
                    2600.0,
                    world.clock.t_myr,
                    6000.0,
                    1.0e-6,
                ) as f32;
                if !d.is_finite() {
                    d = 6000.0;
                }
                world.depth_m[i] = d.clamp(0.0, 6000.0);
            }
            // Build continent template (unitless 0..1)
            let continent_tpl: Vec<f32> = if continents_n == 0 {
                // Supercontinent: contiguous ribbon of lobes
                engine::continent::build_supercontinent_template(
                    &world.grid,
                    ov.simple_seed,
                    5,
                    2200.0,
                    600.0,
                )
            } else {
                let cp = engine::continent::ContinentParams {
                    seed: ov.simple_seed,
                    n_continents: continents_n,
                    mean_radius_km: 2200.0,
                    falloff_km: 600.0,
                    plateau_uplift_m: 1.0,
                    target_land_fraction: None,
                };
                let cf = engine::continent::build_continents(&world.grid, cp);
                cf.uplift_template_m
            };
            // If Supercontinent, imprint broad inherited belts before amplitude solve
            if continents_n == 0 {
                let belts = engine::continent::build_supercontinent_belts(
                    ov.simple_seed,
                    engine::continent::BeltParams {
                        half_width_km_primary: ov.belt_hw_primary_km,
                        half_width_km_secondary: ov.belt_hw_secondary_km,
                        uplift_primary_m: ov.belt_uplift_primary_m,
                        uplift_secondary_m: ov.belt_uplift_secondary_m,
                        diag_angle_deg: ov.belt_diag_deg,
                    },
                );
                let imprint = engine::continent::imprint_orogenic_belts(&world.grid, &belts);
                for (d, di) in world.depth_m.iter_mut().zip(imprint.iter()) {
                    *d += *di;
                }
            }
            // Solve amplitude for target land fraction (area-based). We'll apply uplift only; sea offset handled via world.sea.eta_m
            let target_land = ov.simple_target_land.clamp(0.0, 0.5);
            let (amp_m, _off_m_unused) = engine::isostasy::solve_amplitude_for_land_fraction(
                &continent_tpl,
                &world.depth_m,
                &world.area_m2,
                target_land,
                0.0,
                6000.0,
                2e-2,
                1e-4,
                48,
            );
            // Apply uplift using the solved amplitude; no boost retries
            let amp_used_m = amp_m as f32;
            for (d, t) in world.depth_m.iter_mut().zip(continent_tpl.iter()) {
                *d -= amp_used_m * *t;
            }
            // Solve η so that land fraction with elev=η−depth matches target
            let eta_off = engine::isostasy::solve_offset_for_land_fraction(
                &world.depth_m,
                &world.area_m2,
                ov.simple_target_land,
                64,
            );
            // Solver returns an offset 'off' used in elevation = -(depth + off). Our convention is elev = η − depth, so η = -off.
            world.sea.eta_m = -(eta_off as f32);
            // Seed continents fields so subsequent steps preserve uplift
            world.c.clone_from(&continent_tpl);
            if continents_n == 0 {
                let th = engine::continent::build_craton_thickness(
                    &world.grid,
                    &continent_tpl,
                    22.0,
                    16.0,
                    2.0,
                    ov.simple_seed,
                    1.5,
                );
                world.th_c_m.clone_from_slice(&th);
            } else {
                world.th_c_m.fill(amp_used_m);
            }
            world.epoch_continents = world.epoch_continents.wrapping_add(1);
            // Diagnostics and UI updates using elev = η − depth
            let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
            let mut land_area = 0.0f64;
            for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) {
                if (world.sea.eta_m as f64 - d as f64) > 0.0 {
                    land_area += a as f64;
                }
            }
            let land_frac = if total_area > 0.0 { land_area / total_area } else { 0.0 };
            // No boost retry. Diagnostics only.
            // Elevation stats (m)
            let mut zmin = f32::INFINITY;
            let mut zmax = f32::NEG_INFINITY;
            let mut zsum = 0.0f64;
            let mut zn = 0usize;
            for &d in &world.depth_m {
                if d.is_finite() {
                    let z = world.sea.eta_m - d;
                    zmin = zmin.min(z);
                    zmax = zmax.max(z);
                    zsum += z as f64;
                    zn += 1;
                }
            }
            let _zmean = if zn > 0 { zsum / (zn as f64) } else { 0.0 };
            println!(
                "[simple] target_land={:.2} solved_eta={:+.0} m | land={:.3}",
                ov.simple_target_land, world.sea.eta_m, land_frac
            );
            println!(
                "[draw] elev min/mean/max = {:.0}/{:.0}/{:.0} m | land={:.1}%",
                zmin,
                _zmean,
                zmax,
                land_frac * 100.0
            );
            // Guards (use eta-adjusted elevation for land/ocean checks)
            let has_land = world.depth_m.iter().any(|&d| (world.sea.eta_m as f64 - d as f64) > 0.0);
            let has_ocean =
                world.depth_m.iter().any(|&d| (world.sea.eta_m as f64 - d as f64) <= 0.0);
            debug_assert!(has_land, "no land after solve");
            debug_assert!(has_ocean, "no ocean after solve");
            // Apply palette then mark dirty
            apply_simple_palette(ov, world, ctx);
            ov.world_dirty = true;
            ov.color_dirty = true;
            ov.raster_dirty = true;
            ov.bathy_cache = None;
            ctx.request_repaint();
            // Ready: let the user hit Play to evolve in real-time (non-blocking)
            ov.stepper.playing = false;
            ov.run_active = false;
            ov.run_target_myr = ov.simple_t_end_myr;
            ov.raster_dirty = true;
            ctx.request_repaint();
        }
    });

    egui::CollapsingHeader::new("Continents & Seeds").default_open(true).show(ui, |ui| {
        ui.horizontal(|ui| {
            ui.label("Preset:");
            egui::ComboBox::from_id_source("simple_preset_combo")
                .selected_text(match ov.simple_preset {
                    1 => "Many small",
                    2 => "Few large",
                    3 => "Supercontinent",
                    _ => "Default",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut ov.simple_preset, 0, "Default");
                    ui.selectable_value(&mut ov.simple_preset, 1, "Many small");
                    ui.selectable_value(&mut ov.simple_preset, 2, "Few large");
                    ui.selectable_value(&mut ov.simple_preset, 3, "Supercontinent");
                });
            if ui
                .button("Pangea")
                .on_hover_text("Set Supercontinent preset and regenerate")
                .clicked()
            {
                ov.simple_preset = 3;
                // Force a fresh world regenerate next frame via existing Generate button path
                // by toggling world_dirty and raster_dirty; user clicks Generate World to apply.
                ov.world_dirty = true;
                ov.raster_dirty = true;
            }
        });
        ui.add(egui::DragValue::new(&mut ov.simple_seed).speed(1));
        ui.add(
            egui::Slider::new(&mut ov.simple_target_land, 0.0..=0.6).text("Target land fraction"),
        );
        if ui.button("🎲 Randomise").on_hover_text("Reseed preset RNG").clicked() {
            // reseed preset RNG only
        }
    });

    egui::CollapsingHeader::new("Maps & Colours").default_open(true).show(ui, |ui| {
        let mut changed = false;
        let is_hyp = ov.simple_palette == 0;
        let is_bio = ov.simple_palette == 1;
        if ui.selectable_label(is_hyp, "Hypsometric").clicked() {
            ov.simple_palette = 0;
            changed = true;
        }
        if ui.selectable_label(is_bio, "Biomes").clicked() {
            ov.simple_palette = 1;
            changed = true;
        }
        ui.separator();
        ui.label("GPU raster debug:");
        changed |= ui.checkbox(&mut ov.gpu_dbg_wire, "Wireframe").changed();
        changed |= ui.checkbox(&mut ov.gpu_dbg_face_tint, "Face tint").changed();
        changed |= ui.checkbox(&mut ov.gpu_dbg_grid, "Grid 8x").changed();
        changed |= ui.checkbox(&mut ov.gpu_dbg_tri_parity, "Tri parity").changed();
        changed |= ui.checkbox(&mut ov.gpu_dbg_tri_index, "Tri index").changed();
        // Simple-mode physics isolation toggles (T-646)
        ui.separator();
        ui.checkbox(&mut ov.freeze_eta, "Freeze sea level");
        ui.checkbox(&mut ov.disable_erosion, "Disable erosion");
        ui.checkbox(&mut ov.disable_flexure, "Disable flexure");
        ui.checkbox(&mut ov.disable_subduction, "Disable subduction");
        if changed {
            apply_simple_palette(ov, world, ctx);
        }
        changed |=
            ui.checkbox(&mut ov.dbg_cpu_bary_gpu_lattice, "CPU face -> GPU bary+lattice").changed();
        changed |=
            ui.checkbox(&mut ov.dbg_gpu_face_cpu_lattice, "GPU face+bary -> CPU lattice").changed();
        changed |= ui.checkbox(&mut ov.show_parity_heat, "Parity heat").changed();
        changed |=
            ui.checkbox(&mut ov.force_cpu_face_pick, "Force CPU face pick (debug)").changed();
        changed |=
            ui.checkbox(&mut ov.dbg_rollover_probe, "Log rollover probe CSV (ROI)").changed();
        if changed {
            ov.raster_dirty = true;
        }
        if changed {
            apply_simple_palette(ov, world, ctx);
        }
    });

    egui::CollapsingHeader::new("Export").default_open(true).show(ui, |ui| {
        let export_h = egui::Button::new("Export heightmap PNG");
        let resp_h = ui.add_enabled(false, export_h);
        if resp_h.hovered() {
            resp_h.on_hover_text("Export not wired yet (no image dependency)");
        }
        let export_c = egui::Button::new("Export colour map PNG");
        let resp_c = ui.add_enabled(false, export_c);
        if resp_c.hovered() {
            resp_c.on_hover_text("Export not wired yet (no image dependency)");
        }
        if ui.button("Export parity CSV").clicked() {
            ov.export_parity_csv_requested = true;
        }
    });
}

fn render_advanced_panels(
    ui: &mut egui::Ui,
    _ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
) {
    egui::CollapsingHeader::new("Kinematics & Boundaries")
        .default_open(ov.adv_open_kinematics)
        .show(ui, |ui| {
            let mut changed = false;
            changed |= ui.checkbox(&mut ov.kin_enable, "Enable rigid motion").changed();
            changed |= ui
                .add(egui::Slider::new(&mut ov.kin_trail_steps, 0..=5).text("Trail steps"))
                .changed();
            if changed {
                ov.bounds_cache = None;
                ov.plates_cache = None;
            }
        });

    egui::CollapsingHeader::new("Map Colour").default_open(ov.adv_open_map_color).show(ui, |ui| {
        let mut changed = false;
        let before_mode = ov.color_mode;
        ui.horizontal(|ui| {
            ui.radio_value(&mut ov.color_mode, 0u8, "Hypsometric");
            ui.radio_value(&mut ov.color_mode, 1u8, "Biome preview");
        });
        if ov.color_mode != before_mode {
            changed = true;
        }
        ui.separator();
        changed |= ui.checkbox(&mut ov.shade_on, "Hillshade").changed();
        ui.add_enabled(
            ov.shade_on,
            egui::Slider::new(&mut ov.shade_strength, 0.0..=1.0).text("Strength"),
        );
        ui.horizontal(|ui| {
            ui.add_enabled(
                ov.shade_on,
                egui::Slider::new(&mut ov.sun_az_deg, 0.0..=360.0).text("Sun az (deg)"),
            );
            ui.add_enabled(
                ov.shade_on,
                egui::Slider::new(&mut ov.sun_alt_deg, 0.0..=90.0).text("Sun alt (deg)"),
            );
        });
        changed |= ui.checkbox(&mut ov.legend_on, "Legend").changed();
        if changed {
            ov.color_dirty = true;
            ov.bathy_cache = None;
        }
    });

    egui::CollapsingHeader::new("Flexure").default_open(ov.adv_open_flexure).show(ui, |ui| {
        let mut changed = false;
        changed |= ui.checkbox(&mut ov.enable_flexure, "Enable flexure (apply to depth)").changed();
        changed |= ui.checkbox(&mut ov.show_flexure, "Show w overlay").changed();
        changed |= ui.checkbox(&mut ov.subtract_mean_load, "Subtract mean load").changed();
        changed |= ui.add(egui::Slider::new(&mut ov.e_gpa, 20.0..=120.0).text("E (GPa)")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.nu, 0.15..=0.30).text("nu")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.te_km, 5.0..=50.0).text("Te (km)")).changed();
        changed |=
            ui.add(egui::Slider::new(&mut ov.k_winkler, 0.0..=5.0e8).text("k (N/m^3)")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.wj_omega, 0.6..=0.9).text("ω (WJ)")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.nu1, 0..=4).text("ν1")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.nu2, 0..=4).text("ν2")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.levels, 1..=8).text("Levels")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.flex_cycles, 1..=5).text("V-cycles")).changed();
        changed |= ui
            .add(egui::Slider::new(&mut ov.flex_gain, 1.0..=20.0).text("Overlay gain ×"))
            .changed();
        changed |=
            ui.checkbox(&mut ov.apply_gain_to_depth, "Apply gain to depth (debug)").changed();
        let mut mp = ov.max_points_flex as u32;
        changed |= ui.add(egui::Slider::new(&mut mp, 1000..=50_000).text("Max flex pts")).changed();
        if changed {
            ov.max_points_flex = mp as usize;
            ov.flex_dirty = true;
        }
        ui.label(format!("residual ratio = {:.3}", ov.last_residual));
    });

    egui::CollapsingHeader::new("Surface Processes").default_open(ov.adv_open_surface).show(ui, |ui| {
        ui.checkbox(&mut ov.surface_enable, "Enable erosion & diffusion");
        ui.add(egui::Slider::new(&mut ov.surf_k_stream, 1.0e-7..=1.0e-5).logarithmic(true).text("k_stream"));
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut ov.surf_m_exp, 0.3..=0.8).text("m"));
            ui.add(egui::Slider::new(&mut ov.surf_n_exp, 0.7..=1.5).text("n"));
        });
        ui.add(egui::Slider::new(&mut ov.surf_k_diff, 0.01..=0.5).logarithmic(true).text("κ_diff (m²/yr)"));
        ui.add(egui::Slider::new(&mut ov.surf_k_tr, 0.0..=0.3).logarithmic(true).text("K_transport"));
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut ov.surf_p_exp, 0.8..=1.6).text("p"));
            ui.add(egui::Slider::new(&mut ov.surf_q_exp, 0.5..=1.5).text("q"));
        });
        ui.add(egui::Slider::new(&mut ov.surf_rho_sed, 1200.0..=2600.0).text("ρ_sed (kg/m³)"));
        ui.add(egui::Slider::new(&mut ov.surf_min_slope, 1.0e-5..=1.0e-3).logarithmic(true).text("min slope"));
        let mut sc = ov.surf_subcycles;
        let changed_sc = ui.add(egui::Slider::new(&mut sc, 1..=8).text("Subcycles")).changed();
        if changed_sc { ov.surf_subcycles = sc; }
        ui.add(egui::Slider::new(&mut ov.cadence_surf_every, 1..=20).text("Cadence (steps)"));
        if let Some(stats) = world.last_surface_stats {
            let pct = if stats.eroded_m3 > 0.0 { stats.residual_m3 / stats.eroded_m3 } else { 0.0 };
            ui.label(format!(
                "Erosion={:.2e} m³  Deposition={:.2e} m³  Residual={:+.2}%  max_ero={:.2} m  max_dep={:.2} m",
                stats.eroded_m3, stats.deposited_m3, pct * 100.0, stats.max_erosion_m, stats.max_deposition_m
            ));
        }
    });

    egui::CollapsingHeader::new("Continents").default_open(true).show(ui, |ui| {
        let mut changed = false;
        changed |=
            ui.add(egui::DragValue::new(&mut ov.cont_seed).speed(1.0).prefix("Seed ")).changed();
        changed |= ui.add(egui::Slider::new(&mut ov.cont_n, 1..=6).text("n continents")).changed();
        changed |= ui
            .add(egui::Slider::new(&mut ov.cont_radius_km, 800.0..=3500.0).text("Radius (km)"))
            .changed();
        changed |= ui
            .add(egui::Slider::new(&mut ov.cont_falloff_km, 200.0..=1200.0).text("Falloff σ (km)"))
            .changed();
        ui.separator();
        ui.label("Supercontinent belts (if n=0)");
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.belt_hw_primary_km, 100.0..=800.0)
                    .text("Primary half-width (km)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.belt_hw_secondary_km, 80.0..=600.0)
                    .text("Secondary half-width (km)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.belt_uplift_primary_m, 0.0..=1500.0)
                    .text("Primary uplift (m)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.belt_uplift_secondary_m, 0.0..=1000.0)
                    .text("Secondary uplift (m)"),
            )
            .changed();
        changed |= ui
            .add(egui::Slider::new(&mut ov.belt_diag_deg, 0.0..=60.0).text("Diagonal angle (deg)"))
            .changed();
        changed |= ui.checkbox(&mut ov.cont_auto_amp, "Auto amplitude to target land %").changed();
        if ov.cont_auto_amp {
            changed |= ui
                .add(
                    egui::Slider::new(&mut ov.cont_target_land_frac, 0.10..=0.60)
                        .text("Target land %"),
                )
                .changed();
        } else {
            changed |= ui
                .add(
                    egui::Slider::new(&mut ov.cont_manual_amp_m, 0.0..=5000.0)
                        .text("Manual amplitude (m)"),
                )
                .changed();
        }
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.cont_max_points, 1000..=50_000)
                    .text("Max overlay points"),
            )
            .changed();
        ui.separator();
        ui.label("Plate Type overlay");
        ui.horizontal(|ui| {
            ui.checkbox(&mut ov.show_plate_type, "Show Plate Type");
            let mut mode = ov.plate_type_mode;
            egui::ComboBox::from_id_source("plate_type_mode")
                .selected_text(if mode == 0 { "Plate-level" } else { "Cell-level" })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut mode, 0, "Plate-level");
                    ui.selectable_value(&mut mode, 1, "Cell-level");
                });
            if mode != ov.plate_type_mode {
                ov.plate_type_mode = mode;
            }
        });
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut ov.c_thresh_cont, 0.4..=0.9).text("C_cont_thresh"));
            ui.add(egui::Slider::new(&mut ov.c_thresh_ocean, 0.0..=0.4).text("C_ocean_thresh"));
        });
        ui.label(format!(
            "land = {:.1}%  amplitude = {:.0} m",
            ov.cont_land_frac * 100.0,
            ov.cont_amp_applied_m
        ));
        if changed {
            ov.mesh_continents = None;
            ov.mesh_coastline = None;
            ov.bathy_cache = None;
            ov.world_dirty = true;
        }
    });

    egui::CollapsingHeader::new("Subduction (bands & magnitudes)").default_open(true).show(
        ui,
        |ui| {
            let mut changed = false;
            ui.horizontal(|ui| {
                ui.label("Preset:");
                egui::ComboBox::from_id_source("subduction_preset_combo")
                    .selected_text(match ov.sub_preset {
                        1 => "Strong rollback",
                        2 => "Back-arc extension",
                        3 => "Weak arcs",
                        _ => "Reference",
                    })
                    .show_ui(ui, |ui| {
                        if ui.selectable_label(ov.sub_preset == 0, "Reference").clicked() {
                            ov.sub_preset = 0;
                            apply_sub_preset(ov);
                            changed = true;
                        }
                        if ui.selectable_label(ov.sub_preset == 1, "Strong rollback").clicked() {
                            ov.sub_preset = 1;
                            apply_sub_preset(ov);
                            changed = true;
                        }
                        if ui.selectable_label(ov.sub_preset == 2, "Back-arc extension").clicked() {
                            ov.sub_preset = 2;
                            apply_sub_preset(ov);
                            changed = true;
                        }
                        if ui.selectable_label(ov.sub_preset == 3, "Weak arcs").clicked() {
                            ov.sub_preset = 3;
                            apply_sub_preset(ov);
                            changed = true;
                        }
                    });
            });
            ui.label("Convergence threshold and band geometry");
            changed |= ui
                .add(
                    egui::Slider::new(&mut ov.sub_tau_conv_m_per_yr, 0.001..=0.02)
                        .text("τ_conv (m/yr)"),
                )
                .changed();
            ui.horizontal(|ui| {
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_trench_half_width_km, 10.0..=120.0)
                            .text("Trench half-width (km)"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_arc_offset_km, 60.0..=300.0)
                            .text("Arc offset (km)"),
                    )
                    .changed();
            });
            ui.horizontal(|ui| {
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_arc_half_width_km, 10.0..=80.0)
                            .text("Arc half-width (km)"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_backarc_width_km, 60.0..=300.0)
                            .text("Back-arc width (km)"),
                    )
                    .changed();
            });
            ui.separator();
            ui.label("Magnitudes (m; positive deepens, negative uplifts)");
            changed |= ui
                .add(
                    egui::Slider::new(&mut ov.sub_trench_deepen_m, 500.0..=4000.0)
                        .text("Trench deepen (m)"),
                )
                .changed();
            ui.horizontal(|ui| {
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_arc_uplift_m, -1000.0..=0.0)
                            .text("Arc uplift (m)"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_backarc_uplift_m, -1000.0..=0.0)
                            .text("Back-arc uplift (m)"),
                    )
                    .changed();
            });
            ui.separator();
            ui.label("Rollback and extension mode");
            ui.horizontal(|ui| {
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_rollback_offset_m, 0.0..=200_000.0)
                            .text("Rollback offset (m)"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_rollback_rate_km_per_myr, 0.0..=40.0)
                            .text("Rollback rate (km/Myr)"),
                    )
                    .changed();
            });
            changed |= ui
                .checkbox(&mut ov.sub_backarc_extension_mode, "Back-arc extension mode")
                .changed();
            if ov.sub_backarc_extension_mode {
                changed |= ui
                    .add(
                        egui::Slider::new(&mut ov.sub_backarc_extension_deepen_m, 0.0..=1000.0)
                            .text("Back-arc deepen (m)"),
                    )
                    .changed();
            }
            ui.separator();
            changed |= ui
                .add(
                    egui::Slider::new(&mut ov.sub_continent_c_min, 0.4..=0.9)
                        .text("Continental C threshold"),
                )
                .changed();
            if changed {
                ov.world_dirty = true;
                ov.bathy_cache = None;
            }
        },
    );

    egui::CollapsingHeader::new("Playback & Caps").default_open(ov.adv_open_playback).show(
        ui,
        |ui| {
            ui.horizontal(|ui| {
                ui.add(
                    egui::Slider::new(&mut ov.vel_scale_px_per_cm_yr, 0.1..=2.0)
                        .text("Vel scale (px per cm/yr)"),
                );
                let arrows = egui::Slider::new(&mut ov.max_arrows_slider, 500..=20_000)
                    .text("Max arrows")
                    .step_by(500.0);
                ui.add(arrows);
                let bounds_cap = egui::Slider::new(&mut ov.max_bounds_slider, 500..=20_000)
                    .text("Max boundaries")
                    .step_by(500.0);
                ui.add(bounds_cap);
                let subd_cap = egui::Slider::new(&mut ov.max_subd_slider, 500..=20_000)
                    .text("Max subduction points")
                    .step_by(500.0);
                ui.add(subd_cap);
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.adaptive_cap, "Adaptive cap (16.6 ms target)");
                ui.label(format!(
                    "live arrows={} live boundaries={} live subd={}",
                    ov.live_arrows_cap, ov.live_bounds_cap, ov.live_subd_cap
                ));
            });
        },
    );

    egui::CollapsingHeader::new("Sea Level").default_open(ov.adv_open_sea_level).show(ui, |ui| {
        let mut changed = false;
        changed |= ui.checkbox(&mut ov.apply_sea_level, "Apply global sea level").changed();
        ui.horizontal(|ui| {
            if ui.button("Re-baseline L now").clicked() {
                let area = world.area_m2.clone();
                let _ = engine::isostasy::rebaseline(world, &area);
                ov.world_dirty = true;
                ov.bathy_cache = None;
                ov.raster_dirty = true;
            }
            ui.checkbox(&mut ov.auto_rebaseline_l, "Auto re-baseline after continents change");
        });
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.target_ocean_fraction, 0.05..=0.95)
                    .text("Target ocean fraction")
                    .clamp_to_range(true),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.extra_offset_m, -4000.0..=4000.0)
                    .text("Extra Δoffset (m)"),
            )
            .changed();
        ui.separator();
        ui.checkbox(&mut ov.lock_bathy_scale, "Lock bathy colour scale");
        ui.add(egui::DragValue::new(&mut ov.bathy_min_max.0).speed(10.0).prefix("min "));
        ui.add(egui::DragValue::new(&mut ov.bathy_min_max.1).speed(10.0).prefix("max "));
        if changed {
            ov.bathy_cache = None;
        }
    });

    egui::CollapsingHeader::new("Hypsometry")
 		.default_open(ov.adv_open_hypsometry)
 		.show(ui, |ui| {
			let mut changed = false;
			changed |= ui.add(egui::Slider::new(&mut ov.hyps_bins, 64..=512).text("bins")).changed();
			changed |= ui.checkbox(&mut ov.hyps_auto_domain, "auto domain").changed();
			if !ov.hyps_auto_domain {
				changed |= ui.add(egui::DragValue::new(&mut ov.hyps_min_m).speed(10.0).prefix("min ")).changed();
				changed |= ui.add(egui::DragValue::new(&mut ov.hyps_max_m).speed(10.0).prefix("max ")).changed();
			}
			if ui.button("Export CSV").clicked() {
				let secs = std::time::SystemTime::now()
					.duration_since(std::time::UNIX_EPOCH)
					.map(|d| d.as_secs())
					.unwrap_or(0);
				let name = format!("hypsometry_{}.csv", secs);
				if !ov.hyps_centers_m.is_empty() && ov.hyps_centers_m.len() == ov.hyps_area_per_bin_m2.len() {
					let mut cum: f64 = 0.0;
					let mut s = String::from("bin_center_m,count_cells,area_m2,cumulative_area_m2\n");
					for (i, &c) in ov.hyps_centers_m.iter().enumerate() {
						let a = ov.hyps_area_per_bin_m2[i];
						cum += a;
						s.push_str(&format!("{:.6},{},{} ,{}\n", c, 0, a, cum));
					}
					let _ = std::fs::write(name, s);
				}
			}
			ui.label(format!(
				"[hyps] bins={} land={:.1}% mean={:.0} m median={:.0} m coast_area={:.3e} m^2 | overlay pts={} transforms pts={}",
				ov.hyps_bins, ov.hyps_land_frac * 100.0, ov.hyps_mean_m, ov.hyps_median_m, ov.hyps_coast_area_m2,
				ov.flex_overlay_count, ov.trans_pull_count + ov.trans_rest_count
			));
			if changed { ov.world_dirty = true; }
		});

    egui::CollapsingHeader::new("Transforms").default_open(ov.adv_open_transforms).show(ui, |ui| {
        ui.label("Active if |tangential| ≥ min_tangential & |normal| ≤ τ_open");
        ui.label("Cyan = pull-apart (deeper), Brown = restraining (shallower)");
        ui.label("Width = half-width from the fault trace (km)");
        let mut changed = false;
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.trans_min_tangential_m_per_yr, 0.0001..=0.03)
                    .logarithmic(true)
                    .text("min tangential (m/yr)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.trans_tau_open_m_per_yr, 0.002..=0.02)
                    .logarithmic(true)
                    .text("τ_open (m/yr)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.trans_basin_half_width_km, 10.0..=60.0)
                    .text("Half-width (km)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.trans_basin_deepen_m, 100.0..=1000.0)
                    .text("Basin deepen (m)"),
            )
            .changed();
        changed |= ui
            .add(
                egui::Slider::new(&mut ov.trans_ridge_like_uplift_m, -800.0..=-50.0)
                    .text("Restraining uplift (m)"),
            )
            .changed();
        changed |= ui
            .add(egui::Slider::new(&mut ov.trans_max_points, 500..=20_000).text("Max points"))
            .changed();
        if changed {
            ov.trans_pull = None;
            ov.trans_rest = None;
            ov.bathy_cache = None;
        }
    });

    egui::CollapsingHeader::new("Debug & Cadence").default_open(false).show(ui, |ui| {
        ui.checkbox(&mut ov.debug_enable_all, "Debug: run all passes");
        ui.horizontal(|ui| {
            ui.label("Profiling burst (steps)");
            ui.add(egui::DragValue::new(&mut ov.debug_burst_steps).speed(1));
        });
        ui.separator();
        ui.label("Per-pass cadence (every N steps)");
        ui.horizontal(|ui| {
            ui.label("advection");
            let mut v = ov.cadence_adv_every.max(1);
            if ui.add(egui::DragValue::new(&mut v).clamp_range(1..=u32::MAX).speed(1)).changed() {
                ov.cadence_adv_every = v.max(1);
            }
        });
        ui.horizontal(|ui| {
            ui.label("transforms");
            let mut v = ov.cadence_trf_every.max(1);
            if ui.add(egui::DragValue::new(&mut v).clamp_range(1..=u32::MAX).speed(1)).changed() {
                ov.cadence_trf_every = v.max(1);
            }
        });
        ui.horizontal(|ui| {
            ui.label("subduction");
            let mut v = ov.cadence_sub_every.max(1);
            if ui.add(egui::DragValue::new(&mut v).clamp_range(1..=u32::MAX).speed(1)).changed() {
                ov.cadence_sub_every = v.max(1);
            }
        });
        ui.horizontal(|ui| {
            ui.label("flexure");
            let mut v = ov.cadence_flx_every.max(1);
            if ui.add(egui::DragValue::new(&mut v).clamp_range(1..=u32::MAX).speed(1)).changed() {
                ov.cadence_flx_every = v.max(1);
            }
        });
        ui.horizontal(|ui| {
            ui.label("sea");
            let mut v = ov.cadence_sea_every.max(1);
            if ui.add(egui::DragValue::new(&mut v).clamp_range(1..=u32::MAX).speed(1)).changed() {
                ov.cadence_sea_every = v.max(1);
            }
        });
    });
}

fn apply_simple_palette(
    ov: &mut overlay::OverlayState,
    _world: &mut engine::world::World,
    ctx: &egui::Context,
) {
    ov.color_mode = if ov.simple_palette == 0 { 0 } else { 1 };
    ov.show_bathy = true;
    ov.bathy_cache = None;
    ov.color_dirty = true;
    ctx.request_repaint();
}
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
    depth_tex: wgpu::Texture,
    depth_view: wgpu::TextureView,
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

        let depth_tex = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("depth tex"),
            size: wgpu::Extent3d {
                width: config.width,
                height: config.height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_tex.create_view(&wgpu::TextureViewDescriptor::default());

        Self { _instance: instance, surface, device, queue, config, depth_tex, depth_view }
    }

    fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
            self.depth_tex = self.device.create_texture(&wgpu::TextureDescriptor {
                label: Some("depth tex"),
                size: wgpu::Extent3d {
                    width: self.config.width,
                    height: self.config.height,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: wgpu::TextureFormat::Depth32Float,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
                view_formats: &[],
            });
            self.depth_view = self.depth_tex.create_view(&wgpu::TextureViewDescriptor::default());
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
    // Compute raster shader (WGSL)
    let raster_shader = gpu.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("raster shader"),
        source: wgpu::ShaderSource::Wgsl(include_str!("../../shaders/raster.wgsl").into()),
    });
    struct FaceCache {
        f: u32,
        #[allow(dead_code)]
        face_ids: Vec<u32>,
        #[allow(dead_code)]
        face_offs: Vec<u32>,
        #[allow(dead_code)]
        face_geom: Vec<[f32; 4]>,
    }
    struct AppPipelines {
        raster: Option<raster_gpu::RasterGpu>,
        face_cache: Option<FaceCache>,
        globe_mesh: Option<globe::GlobeMesh>,
        globe: Option<globe::GlobeRenderer>,
        globe_cam: Option<globe::OrbitCamera>,
    }
    let mut pipes = AppPipelines {
        raster: None,
        face_cache: None,
        globe_mesh: None,
        globe: None,
        globe_cam: None,
    };
    let mut ov = overlay::OverlayState::default();
    // Edge-triggered snapshots state
    let mut next_snapshot_t: f64 = f64::INFINITY;
    let mut flex = plot_flexure::FlexureUI::default();
    let mut age_plot = plot_age_depth::AgeDepthUIState::default();
    // T-020: Construct device field buffers sized to the grid (then drop)
    {
        let f: u32 = 64;
        let g_tmp = engine::grid::Grid::new(f);
        let _device_fields = engine::fields::DeviceFields::new(&gpu.device, g_tmp.cells);
    }
    // Obsolete: replaced by ov.view_mode
    let _globe_enabled: bool = true;
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
    let _nplates: usize = world.plates.pole_axis.len();

    // Playback state
    let playing: bool = false;
    let dt_myr: f32 = 1.0;
    let steps_per_sec: u32 = 5;
    let mut step_accum_s: f32 = 0.0;
    // Snapshots frequency (Myr)
    let snapshot_every_myr: f32 = 5.0;

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
                                let mut _continents_dirty: bool = false;
                                let mut _flex_dirty: bool = false;
                                if !ov.t505_logged {
                                    println!("[ui] T-505 done | drawer={} | mode={}", ov.drawer_open, if ov.mode_simple { "simple" } else { "advanced" });
                                    ov.t505_logged = true;
                                }
                            if !ov.mode_simple {
                                if ctx.input(|i| i.key_pressed(egui::Key::Num1)) { ov.show_plates = !ov.show_plates; ov.plates_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num2)) { ov.show_vel = !ov.show_vel; ov.vel_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num3)) { ov.show_bounds = !ov.show_bounds; ov.bounds_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num4)) { ov.show_plate_adjacency = !ov.show_plate_adjacency; ov.net_adj_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num5)) { ov.show_triple_junctions = !ov.show_triple_junctions; ov.net_tj_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num4)) { ov.show_age = !ov.show_age; ov.age_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num5)) { ov.show_bathy = !ov.show_bathy; ov.bathy_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num6)) { ov.show_age_depth = !ov.show_age_depth; }
                                if ctx.input(|i| i.key_pressed(egui::Key::A)) { age_plot.show = !age_plot.show; }
                                if ctx.input(|i| i.key_pressed(egui::Key::M)) { ov.show_map_color_panel = !ov.show_map_color_panel; }
                                if ctx.input(|i| i.key_pressed(egui::Key::S)) { ov.surface_enable = !ov.surface_enable; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num7)) { ov.show_subduction = !ov.show_subduction; ov.subd_trench=None; ov.subd_arc=None; ov.subd_backarc=None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num0)) { ov.show_transforms = !ov.show_transforms; ov.trans_pull=None; ov.trans_rest=None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::C)) { ov.show_continents = !ov.show_continents; if ov.show_continents && (ov.mesh_continents.is_none() || ov.mesh_coastline.is_none()) { _continents_dirty = true; } }
                                if ctx.input(|i| i.key_pressed(egui::Key::L)) { ov.apply_sea_level = !ov.apply_sea_level; ov.bathy_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::F)) { flex.show = !flex.show; if flex.show { flex.recompute(); } }
                                if ctx.input(|i| i.key_pressed(egui::Key::G)) { ov.show_flexure = !ov.show_flexure; _flex_dirty = true; }
                                if ctx.input(|i| i.key_pressed(egui::Key::H)) { ov.show_hud = !ov.show_hud; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Y)) { ov.show_hypsometry = !ov.show_hypsometry; }
                            } else {
                                // In Simple mode, force color layer on and ignore hide-map hotkeys
                                ov.show_bathy = true;
                            }

                            egui::TopBottomPanel::top("hud").show(ctx, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(format!("Aulé Viewer v{}", env!("CARGO_PKG_VERSION")));
                                    ui.separator();
                                    ui.label("Mode:");
                                    ui.checkbox(&mut ov.mode_simple, "Simple mode");
                                    ui.separator();
                                    ui.label("View:");
                                    ui.selectable_value(&mut ov.view_mode, overlay::ViewMode::Map, "Map 2D");
                                    ui.selectable_value(&mut ov.view_mode, overlay::ViewMode::Globe, "Globe 3D");
                                    if ov.view_mode == overlay::ViewMode::Globe {
                                        ui.add(egui::Slider::new(&mut ov.globe.exaggeration, 0.0..=5.0).text("Z exaggeration"));
                                        ui.checkbox(&mut ov.globe.show_wireframe, "Wireframe");
                                        ui.checkbox(&mut ov.globe.show_probe, "Hover probe");
                                    }
                                    if ui.button(if ov.drawer_open { "⟨⟩" } else { "☰" }).clicked() {
                                        ov.drawer_open = !ov.drawer_open;
                                    }
                                    ui.separator();
                                    // Live land fraction (area-weighted) using elev = η − depth
                                    let area_total: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                                    let mut land_area: f64 = 0.0;
                                    for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) {
                                        if (world.sea.eta_m as f64 - d as f64) > 0.0 { land_area += a as f64; }
                                    }
                                    let land_frac = if area_total > 0.0 { land_area / area_total } else { 0.0 };
                                    ui.label(format!("Land: {:.1}%", land_frac * 100.0));
                                });
                            });
                            // Left drawer under the top bar
                            if ov.drawer_open {
                                egui::SidePanel::left("drawer")
                                .resizable(true)
                                .min_width(280.0)
                                .default_width(320.0)
                                .show(ctx, |ui| {
                                    egui::ScrollArea::vertical()
                                        .auto_shrink([false, false])
                                        .show(ui, |ui| {
                                                if ov.mode_simple {
                                                    render_simple_panels(ui, ctx, &mut world, &mut ov);
                                                } else {
                                                    render_advanced_panels(ui, ctx, &mut world, &mut ov);
                                                }
                                                ui.separator();
                                                ui.collapsing("Debug", |ui| {
                                                    if ui.button("Export raster debug CSV").clicked() {
                                                        let w = 256u32;
                                                        let h = 128u32;
                                                        if let Err(e) = viewer::export::export_raster_debug_csv(&world, w, h, "raster_debug.csv") {
                                                            eprintln!("[debug] export_raster_debug_csv failed: {}", e);
                                                        } else {
                                                            println!("[debug] wrote raster_debug.csv ({}x{})", w, h);
                                                        }
                                                    }
                                                });
                                        });
                                });
                            }
                            // Central canvas (draw only, no controls) — make transparent so 3D pass remains visible
                            egui::CentralPanel::default().frame(egui::Frame::none()).show(ctx, |ui| {
                                let rect = ui.max_rect();
                                let painter = ui.painter_at(rect);
                                if ov.view_mode == overlay::ViewMode::Globe {
                                    // Build globe once (or when F changes later)
                                    if pipes.globe_mesh.is_none() || pipes.globe.is_none() {
                                        let f_now = world.grid.freq();
                                        let mesh = globe::build_globe_mesh(&gpu.device, &world.grid);
                                        let gr = globe::GlobeRenderer::new(&gpu.device, gpu.config.format, mesh.vertex_count);
                                        // Upload initial heights and LUT
                                        gr.upload_heights(&gpu.queue, &world.depth_m);
                                        gr.write_lut_from_overlay(&gpu.queue, &ov);
                                        pipes.globe_mesh = Some(mesh);
                                        pipes.globe = Some(gr);
                                        pipes.globe_cam = Some(globe::OrbitCamera::default());
                                        let _ = f_now; // silence if unused in some cfgs
                                    }

                                    // Update camera from input
                                    if let Some(cam) = &mut pipes.globe_cam {
                                        cam.aspect = (gpu.config.width.max(1) as f32) / (gpu.config.height.max(1) as f32);
                                        let ui_hijacked = ctx.is_using_pointer() || ctx.is_pointer_over_area();
                                        cam.update_from_input(&egui_ctx, ui_hijacked);
                                    }
                                } else if ov.mode_simple {
                                    // T-902A-GPU: compute raster path
                                    if ov.use_gpu_raster {
                                        // Resolution policy: HQ when paused, LQ when running
                                        if !ov.run_active && ov.high_quality_when_paused {
                                            if ov.raster_size != (2048, 1024) { ov.raster_size = (2048, 1024); ov.raster_dirty = true; }
                                        } else if ov.raster_size != (1024, 512) { ov.raster_size = (1024, 512); ov.raster_dirty = true; }
                                        let (rw, rh) = ov.raster_size;
                                        let need_new = pipes
                                            .raster
                                            .as_ref()
                                            .map(|r| r.width != rw || r.height != rh)
                                            .unwrap_or(true);
                                        if need_new {
                                            let rg = raster_gpu::RasterGpu::new(&gpu.device, &raster_shader, rw, rh);
                                            pipes.raster = Some(rg);
                                            ov.raster_dirty = true;
                                            ov.raster_tex_id = None;
                                        }
                                        if let Some(rg) = pipes.raster.as_mut() {
                                            let f_now = world.grid.frequency;
                                            // Build/cache face tables once per F
                                            let cache_stale = match &pipes.face_cache { Some(fc) => fc.f != f_now, None => true };
                                            if cache_stale || rg.face_count == 0 {
                                                let (face_ids_ref, face_offs_ref) = world.grid.face_vertex_table();
                                                let face_ids: Vec<u32> = face_ids_ref.to_vec();
                                                let face_offs: Vec<u32> = face_offs_ref.to_vec();
                                                // Build canonical face geometry, neighbors, and corner permutations from shared crate
                                                let faces = aule_geo::build_face_table();
                                                let (gpu_faces, gpu_neighbors) = aule_geo::to_gpu_faces(&faces);
                                                let perms = aule_geo::corner_permutation_opp(&faces);
                                                let gpu_perms = aule_geo::to_gpu_perms(&perms);
                                                let mut face_geom: Vec<[f32; 4]> = Vec::with_capacity(20 * 4);
                                                for gf in &gpu_faces {
                                                    face_geom.push(gf.a);
                                                    face_geom.push(gf.b);
                                                    face_geom.push(gf.c);
                                                    face_geom.push(gf.n);
                                                }
                                                // Build neighbor mapping from shared crate (opp A,B,C). Shader uses neighbor ids and perms.
                                                let mut face_edge_info: Vec<u32> = Vec::with_capacity(20 * 3 * 3);
                                                for n in &gpu_neighbors {
                                                    face_edge_info.push(n.opp[0]); face_edge_info.push(0); face_edge_info.push(0);
                                                    face_edge_info.push(n.opp[1]); face_edge_info.push(0); face_edge_info.push(0);
                                                    face_edge_info.push(n.opp[2]); face_edge_info.push(0); face_edge_info.push(0);
                                                }
                                                // Append perms as tightly packed u32 triplets (3 edges × 3 entries), each entry packs 3 u8: [permA,permB,permC]
                                                let mut face_perm_info: Vec<u32> = Vec::with_capacity(20 * 3);
                                                for gp in &gpu_perms {
                                                    for e in 0..3 {
                                                        let p = gp.perm_opp[e];
                                                        let packed = (p[0] as u32) | ((p[1] as u32) << 8) | ((p[2] as u32) << 16);
                                                        face_perm_info.push(packed);
                                                    }
                                                }
                                                let dbg = (if ov.gpu_dbg_wire { 1u32 } else { 0 })
                                                    | (if ov.gpu_dbg_face_tint { 1u32<<1 } else { 0 })
                                                    | (if ov.gpu_dbg_grid { 1u32<<2 } else { 0 })
                                                    | (if ov.gpu_dbg_tri_parity { 1u32<<3 } else { 0 })
                                                    | (if ov.gpu_dbg_tri_index { 1u32<<5 } else { 0 })
                                                    | (if ov.show_parity_heat || ov.export_parity_csv_requested { 1u32<<6 } else { 0 })
                                                    | (if ov.force_cpu_face_pick { 1u32<<7 } else { 0 })
                                                    | (if ov.dbg_cpu_bary_gpu_lattice { 1u32<<10 } else { 0 });
                                                let u = raster_gpu::Uniforms { width: rw, height: rh, f: f_now, palette_mode: if ov.color_mode == 0 { 0 } else { 1 }, debug_flags: dbg, d_max: ov.hypso_d_max.max(1.0), h_max: ov.hypso_h_max.max(1.0), snowline: ov.hypso_snowline, eta_m: world.sea.eta_m, inv_dmax: 1.0f32/ov.hypso_d_max.max(1.0), inv_hmax: 1.0f32/ov.hypso_h_max.max(1.0) };
                                                rg.upload_inputs(&gpu.device, &gpu.queue, &u, face_ids.clone(), face_offs.clone(), &face_geom, &world.depth_m, &face_edge_info, &face_perm_info);
                                                rg.write_lut_from_overlay(&gpu.queue, &ov);
                                                // If forcing CPU face pick, generate per-pixel face ids using shared picker
                                                if ov.force_cpu_face_pick || ov.dbg_cpu_bary_gpu_lattice {
                                                    let mut cpu_faces: Vec<u32> = vec![0; (rw * rh) as usize];
                                                    let picker = GeoPicker::new();
                                                    for y in 0..rh { for x in 0..rw {
                                                        let (lon, lat) = pixel_to_lon_lat(x, y, rw, rh);
                                                        let p3 = sph_to_unit_px(lon, lat);
                                                        let p = Vec3::new(p3.x, p3.y, p3.z).norm();
                                                        let res = picker.pick_from_unit(p, f_now);
                                                        cpu_faces[(y*rw + x) as usize] = res.face;
                                                    }}
                                                    rg.write_cpu_face_pick(&gpu.queue, &cpu_faces);
                                                }
                                                pipes.face_cache = Some(FaceCache { f: f_now, face_ids, face_offs, face_geom });
                                                // Register texture once
                                                if ov.raster_tex_id.is_none() {
                                                    let tid = egui_renderer.register_native_texture(&gpu.device, &rg.out_view, wgpu::FilterMode::Linear);
                                                    ov.raster_tex_id = Some(tid);
                                                }
                                                ov.raster_dirty = false; ov.world_dirty = false; ov.color_dirty = false; ov.last_raster_at = std::time::Instant::now();
                                                println!("[viewer] raster(gpu) W={} H={} | F={} | verts={} | face_tbl={} | dispatch={}x{}", rw, rh, f_now, world.depth_m.len(), 20 * ((f_now + 1) * (f_now + 2) / 2), (rw + 7) / 8, (rh + 7) / 8);
                                            } else {
                                                // Frame-driven stepping (non-blocking)
                                                if ov.stepper.playing && (world.clock.t_myr as f32) < ov.stepper.t_target_myr {
                                                    let start = std::time::Instant::now();
                                                let max_n = ov.stepper.max_steps_per_frame.max(1);
                                                for _ in 0..max_n {
                                                        if (world.clock.t_myr as f32) >= ov.stepper.t_target_myr { break; }
                                                        let _burst = ov.debug_burst_steps > 0;
                                                        let _enable_all = ov.debug_enable_all || _burst;
                                                        let _adv_every = ov.cadence_adv_every.max(1);
                                                        let _trf_every = ov.cadence_trf_every.max(1);
                                                        let _sub_every = ov.cadence_sub_every.max(1);
                                                        let _flx_every = ov.cadence_flx_every.max(1);
                                                        let _sea_every = ov.cadence_sea_every.max(1);
                                                        // (legacy StepParams block removed; pipeline used below)
                                                        let t0 = world.clock.t_myr;
                                                        let cfg = engine::pipeline::PipelineCfg {
                                                            dt_myr,
                                                            steps_per_frame: 1,
                                                            enable_flexure: !ov.disable_flexure,
                                                            enable_erosion: !ov.disable_erosion,
                                                            target_land_frac: ov.simple_target_land,
                                                            freeze_eta: ov.freeze_eta,
                                                            log_mass_budget: false,
                                                            enable_subduction: !ov.disable_subduction,
                                                            enable_rigid_motion: true,
                                                            cadence_trf_every: ov.cadence_trf_every.max(1),
                                                            cadence_sub_every: ov.cadence_sub_every.max(1),
                                                            cadence_flx_every: ov.cadence_flx_every.max(1),
                                                            cadence_sea_every: ov.cadence_sea_every.max(1),
                                                            cadence_surf_every: ov.cadence_sea_every.max(1),
                                                            substeps_transforms: 4,
                                                            substeps_subduction: 4,
                                                            use_gpu_flexure: false,
                                                            gpu_flex_levels: ov.levels.max(1),
                                                            gpu_flex_cycles: ov.flex_cycles.max(1),
                                                            gpu_wj_omega: ov.wj_omega,
                                                            subtract_mean_load: ov.subtract_mean_load,
                                                            surf_k_stream: ov.surf_k_stream,
                                                            surf_m_exp: ov.surf_m_exp,
                                                            surf_n_exp: ov.surf_n_exp,
                                                            surf_k_diff: ov.surf_k_diff,
                                                            surf_k_tr: ov.surf_k_tr,
                                                            surf_p_exp: ov.surf_p_exp,
                                                            surf_q_exp: ov.surf_q_exp,
                                                            surf_rho_sed: ov.surf_rho_sed,
                                                            surf_min_slope: ov.surf_min_slope,
                                                            surf_subcycles: ov.surf_subcycles.max(1),
                                                            surf_couple_flexure: ov.surf_couple_flexure,
                                                            sub_tau_conv_m_per_yr: ov.sub_tau_conv_m_per_yr,
                                                            sub_trench_half_width_km: ov.sub_trench_half_width_km,
                                                            sub_arc_offset_km: ov.sub_arc_offset_km,
                                                            sub_arc_half_width_km: ov.sub_arc_half_width_km,
                                                            sub_backarc_width_km: ov.sub_backarc_width_km,
                                                            sub_trench_deepen_m: ov.sub_trench_deepen_m,
                                                            sub_arc_uplift_m: ov.sub_arc_uplift_m,
                                                            sub_backarc_uplift_m: ov.sub_backarc_uplift_m,
                                                            sub_rollback_offset_m: ov.sub_rollback_offset_m,
                                                            sub_rollback_rate_km_per_myr: ov.sub_rollback_rate_km_per_myr,
                                                            sub_backarc_extension_mode: ov.sub_backarc_extension_mode,
                                                            sub_backarc_extension_deepen_m: ov.sub_backarc_extension_deepen_m,
                                                            sub_continent_c_min: ov.sub_continent_c_min,
                                                            cadence_spawn_plate_every: 0,
                                                            cadence_retire_plate_every: 0,
                                                            cadence_force_balance_every: 8,
                                                            fb_gain: 1.0e-12,
                                                            fb_damp_per_myr: 0.2,
                                                            fb_k_conv: 1.0,
                                                            fb_k_div: 0.5,
                                                            fb_k_trans: 0.1,
                                                            fb_max_domega: 5.0e-9,
                                                            fb_max_omega: 2.0e-7,
                                                        };
                                                        let mut elev_tmp: Vec<f32> = vec![0.0; world.depth_m.len()];
                                                        let mut eta = world.sea.eta_m;
                                                        engine::pipeline::step_full(
                                                            &mut world,
                                                            engine::pipeline::SurfaceFields { elev_m: &mut elev_tmp, eta_m: &mut eta },
                                                            cfg,
                                                        );
                                                        world.sea.eta_m = eta;
                                                        // Elevation stats and logging: z in elev_tmp is already (η − depth); do not add η again
                                                        // Invalidate plate/boundary caches so overlays reflect updated plate_id and edges
                                                        ov.plates_cache = None;
                                                        ov.bounds_cache = None;
                                                        let mut zmin = f32::INFINITY;
                                                        let mut zmax = f32::NEG_INFINITY;
                                                        let mut zsum = 0.0f64;
                                                        let mut land_area = 0.0f64;
                                                        let mut total_area = 0.0f64;
                                                        for (&z0, &a) in elev_tmp.iter().zip(world.area_m2.iter()) {
                                                             if z0.is_finite() {
                                                                 let z = z0; // already η − depth
                                                                 zmin = zmin.min(z);
                                                                 zmax = zmax.max(z);
                                                                 zsum += z as f64;
                                                                 total_area += a as f64;
                                                                 if z > 0.0 { land_area += a as f64; }
                                                             }
                                                         }
                                                        let zmean = if total_area > 0.0 { zsum / total_area } else { 0.0 };
                                                        println!("[draw] elev min/mean/max = {:.0}/{:.0}/{:.0} m | land={:.1}%", zmin, zmean, zmax, 100.0 * land_area / total_area.max(1.0));
                                                        // Bookkeeping
                                                        ov.color_dirty = true;
                                                        ov.world_dirty = true;
                                                        ctx.request_repaint();
                                                        // Per-step land fraction (area-weighted) using elev = η − depth
                                                        let area_total: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                                                        let mut land_area: f64 = 0.0;
                                                        for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) {
                                                            if (world.sea.eta_m as f64 - d as f64) > 0.0 { land_area += a as f64; }
                                                        }
                                                        let land_frac = if area_total > 0.0 { land_area / area_total } else { 0.0 };
                                                        println!(
                                                            "[step/ui] t={:.1}→{:.1} Myr | land={:.1}%",
                                                            t0,
                                                            world.clock.t_myr,
                                                            land_frac * 100.0
                                                        );
                                                        ov.world_dirty = true;
                                                        if start.elapsed().as_millis() > 10 { break; }
                                                    }
                                                    if ov.debug_burst_steps > 0 { ov.debug_burst_steps = ov.debug_burst_steps.saturating_sub(max_n); if ov.debug_burst_steps == 0 { println!("[debug] profiling burst complete."); } }
                                                    // GPU raster throttle while playing
                                                    let now = std::time::Instant::now();
                                                    if now.duration_since(ov.stepper.last_raster_at).as_millis() as u64 >= ov.stepper.min_ms_between_raster {
                                                        ov.raster_dirty = true;
                                                        ov.stepper.last_raster_at = now;
                                                    }
                                                    ctx.request_repaint();
                                                }
                                                // Dispatch only when flagged dirty
                                                if ov.raster_dirty {
                                                    let mut dbg = if ov.gpu_dbg_wire { 1u32 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_face_tint { 1u32<<1 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_grid { 1u32<<2 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_tri_parity { 1u32<<3 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_tri_index { 1u32<<5 } else { 0 };
                                                    dbg |= if ov.show_parity_heat || ov.export_parity_csv_requested { 1u32<<6 } else { 0 };
                                                    dbg |= if ov.force_cpu_face_pick { 1u32<<7 } else { 0 };
                                                    dbg |= if ov.dbg_cpu_bary_gpu_lattice { 1u32<<10 } else { 0 };
                                                    let u = raster_gpu::Uniforms { width: rw, height: rh, f: f_now, palette_mode: if ov.color_mode == 0 { 0 } else { 1 }, debug_flags: dbg, d_max: ov.hypso_d_max.max(1.0), h_max: ov.hypso_h_max.max(1.0), snowline: ov.hypso_snowline, eta_m: world.sea.eta_m, inv_dmax: 1.0f32/ov.hypso_d_max.max(1.0), inv_hmax: 1.0f32/ov.hypso_h_max.max(1.0) };
                                                    // If forcing CPU face pick, generate per-pixel face ids using FACE_GEOM (A,B,C,N)
                                                    if ov.force_cpu_face_pick || ov.dbg_cpu_bary_gpu_lattice {
                                                        let mut cpu_faces: Vec<u32> = vec![0; (rw * rh) as usize];
                                                        let picker = GeoPicker::new();
                                                        for y in 0..rh { for x in 0..rw {
                                                            let (lon, lat) = pixel_to_lon_lat(x, y, rw, rh);
                                                            let p3 = sph_to_unit_px(lon, lat);
                                                            let p = Vec3::new(p3.x, p3.y, p3.z).norm();
                                                            let res = picker.pick_from_unit(p, f_now);
                                                            cpu_faces[(y*rw + x) as usize] = res.face;
                                                        }}
                                                        rg.write_cpu_face_pick(&gpu.queue, &cpu_faces);
                                                    }
                                                    rg.write_uniforms(&gpu.queue, &u);
                                                    rg.write_vertex_values(&gpu.queue, &world.depth_m);
                                                    rg.write_lut_from_overlay(&gpu.queue, &ov);
                                                    rg.dispatch(&gpu.device, &gpu.queue);
                                                    // Optional parity readback and console stats when debug bit is set
                                                    if (dbg & (1u32<<6)) != 0 {
                                                        if let Some(dbg_buf) = rg.read_debug_face_tri(&gpu.device, &gpu.queue) {
                                                            // Sample every 4th pixel with the shared GeoPicker to mirror WGSL exactly
                                                            let stride = 4u32;
                                                            let mut mismatches: u64 = 0;
                                                            let mut face_mismatches: u64 = 0;
                                                            let mut tri_mismatches_same_face: u64 = 0;
                                                            let mut total: u64 = 0;
                                                            let mut pts: Vec<[u32;2]> = Vec::new();
                                                            let (w, h, f) = (rw, rh, f_now);
                                                            let picker = GeoPicker::new();
                                                            for y in (0..h).step_by(stride as usize) {
                                                                for x in (0..w).step_by(stride as usize) {
                                                                    let idx = ((y * w + x) * 2) as usize;
                                                                    let face_gpu = dbg_buf[idx];
                                                                    let tri_gpu_global = dbg_buf[idx + 1];
                                                                    let fx = x as f32 + 0.5; let fy = y as f32 + 0.5;
                                                                    let lon = (fx / w as f32) * std::f32::consts::TAU - std::f32::consts::PI;
                                                                    let lat = std::f32::consts::FRAC_PI_2 - (fy / h as f32) * std::f32::consts::PI;
                                                                    let cl = lat.cos(); let p = Vec3::new(cl*lon.cos(), lat.sin(), -cl*lon.sin()).norm();
                                                                    let cpu = picker.pick_from_unit(p, f);
                                                                    let tri_cpu_global = cpu.face * f * f + cpu.tri;
                                                                    if face_gpu != cpu.face { face_mismatches += 1; mismatches += 1; }
                                                                    else if tri_gpu_global != tri_cpu_global { tri_mismatches_same_face += 1; mismatches += 1; if ov.show_parity_heat { pts.push([x, y]); } }
                                                                    total += 1;
                                                                }
                                                            }
                                                            let pct = if total>0 { 100.0*(mismatches as f64)/(total as f64) } else { 0.0 };
                                                            let pct_face = if total>0 { 100.0*(face_mismatches as f64)/(total as f64) } else { 0.0 };
                                                            let pct_tri = if total>0 { 100.0*(tri_mismatches_same_face as f64)/(total as f64) } else { 0.0 };
                                                            println!("[parity] sampled={} mismatch% total={:.3} face={:.3} tri(sameface)={:.3}", total, pct, pct_face, pct_tri);
                                                            if ov.show_parity_heat { ov.parity_points = Some(pts); }
                                                            // ROI diagnostics: dump a small window with raw barycentrics and neighbor choice (overwrite each dispatch)
                                                            if let Some(fc) = pipes.face_cache.as_ref() {
                                                                use std::io::Write;
                                                                if let Ok(mut fcsv) = std::fs::File::create("parity_roi.csv") {
                                                                    let _ = writeln!(fcsv, "x,y,face_gpu,face_sim,face_cpu,kneg_sim,kneg_cpu,wa_s,wb_s,wc_s,wa_c,wb_c,wc_c");
                                                                    // Build neighbor table once
                                                                    let mut neighbors: [[u32;3]; 20] = [[u32::MAX;3]; 20];
                                                                    let corners = world.grid.face_corners();
                                                                    for fid in 0..20usize {
                                                                        let tri = corners[fid];
                                                                        let mk = |va:u32, vb:u32| -> u32 {
                                                                            let mut nf=u32::MAX;
                                                                            's: for (j,t) in corners.iter().enumerate(){
                                                                                if j==fid {continue;}
                                                                                let arr=*t;
                                                                                let mut a0=None;
                                                                                let mut a1=None;
                                                                                for k in 0..3u32 {
                                                                                    if arr[k as usize]==va { a0=Some(k); }
                                                                                    if arr[k as usize]==vb { a1=Some(k); }
                                                                                }
                                                                                if a0.is_some() && a1.is_some(){ nf=j as u32; break 's; }
                                                                            }
                                                                            nf
                                                                        };
                                                                        neighbors[fid][0] = mk(tri[1], tri[2]);
                                                                        neighbors[fid][1] = mk(tri[2], tri[0]);
                                                                        neighbors[fid][2] = mk(tri[0], tri[1]);
                                                                    }
                                                                    let y0 = (h/3).max(1); let x0 = (w/3).max(1);
                                                                    let hroi = 64u32.min(h); let wroi = 128u32.min(w);
                                                                    for dy in 0..hroi { for dx in 0..wroi {
                                                                        let x = x0 + dx; let y = y0 + dy; let idx = ((y * w + x) * 2) as usize;
                                                                        let face_gpu = dbg_buf[idx];
                                                                        let fx = x as f32 + 0.5; let fy = y as f32 + 0.5;
                                                                        let lon = (fx / w as f32) * std::f32::consts::TAU - std::f32::consts::PI;
                                                                        let lat = std::f32::consts::FRAC_PI_2 - (fy / h as f32) * std::f32::consts::PI;
                                                                        let cl = lat.cos(); let p = [cl*lon.cos(), lat.sin(), -cl*lon.sin()];
                                                                        // Sim path using FACE_GEOM
                                                                        let mut best_f = 0usize; let mut bd = -1e9f32;
                                                                        for fidx in 0..20usize { let n = fc.face_geom[fidx*4+3]; let d = n[0]*p[0]+n[1]*p[1]+n[2]*p[2]; if d>bd { bd=d; best_f=fidx; } }
                                                                        let mut f_sim = best_f as u32; let geom = &fc.face_geom[(f_sim as usize)*4 .. (f_sim as usize)*4 + 4];
                                                                        let a = geom[0]; let b = geom[1]; let c = geom[2]; let n = geom[3];
                                                                        let t = n[0]*(p[0]-a[0]) + n[1]*(p[1]-a[1]) + n[2]*(p[2]-a[2]);
                                                                        let q = [p[0]-n[0]*t, p[1]-n[1]*t, p[2]-n[2]*t];
                                                                        let v0 = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]; let v1 = [c[0]-a[0], c[1]-a[1], c[2]-a[2]]; let v2 = [q[0]-a[0], q[1]-a[1], q[2]-a[2]];
                                                                        let d00 = v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]; let d01 = v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2]; let d11 = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]; let d20 = v2[0]*v0[0]+v2[1]*v0[1]+v2[2]*v0[2]; let d21 = v2[0]*v1[0]+v2[1]*v1[1]+v2[2]*v1[2];
                                                                        let denom = (d00*d11 - d01*d01).max(1e-12);
                                                                        let wb = (d11*d20 - d01*d21)/denom; let wc = (d00*d21 - d01*d20)/denom; let wa = 1.0 - wb - wc;
                                                                        let mut kneg_s = 0usize; let mut vmin = wa;
                                                                        if wb < vmin { vmin = wb; kneg_s = 1; }
                                                                        if wc < vmin { kneg_s = 2; }
                                                                        if wa < -1e-6 || wb < -1e-6 || wc < -1e-6 { let nf = neighbors[f_sim as usize][kneg_s]; if nf!=u32::MAX { f_sim = nf; } }
                                                                        // CPU-grid path
                                                                        let mut best_fg = 0usize; let mut bdg = -1e9f32;
                                                                        for (fi, tri) in corners.iter().enumerate().take(20) {
                                                                            let a = world.grid.pos_xyz[tri[0] as usize]; let b = world.grid.pos_xyz[tri[1] as usize]; let c = world.grid.pos_xyz[tri[2] as usize];
                                                                            let n = [ (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]), (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]), (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) ];
                                                                            let d = n[0]*p[0] + n[1]*p[1] + n[2]*p[2]; if d>bdg { bdg=d; best_fg=fi; }
                                                                        }
                                                                        let mut f_cpu = best_fg as u32; let tri = corners[f_cpu as usize];
                                                                        let a2 = world.grid.pos_xyz[tri[0] as usize]; let b2 = world.grid.pos_xyz[tri[1] as usize]; let c2 = world.grid.pos_xyz[tri[2] as usize];
                                                                        let ab = [b2[0]-a2[0],b2[1]-a2[1],b2[2]-a2[2]]; let ac = [c2[0]-a2[0],c2[1]-a2[1],c2[2]-a2[2]];
                                                                        let nx = ab[1]*ac[2]-ab[2]*ac[1]; let ny = ab[2]*ac[0]-ab[0]*ac[2]; let nz = ab[0]*ac[1]-ab[1]*ac[0]; let inv = 1.0f32/(nx*nx+ny*ny+nz*nz).sqrt().max(1e-8); let n2 = [nx*inv, ny*inv, nz*inv];
                                                                        let t2 = n2[0]*(p[0]-a2[0])+n2[1]*(p[1]-a2[1])+n2[2]*(p[2]-a2[2]); let q2 = [p[0]-n2[0]*t2,p[1]-n2[1]*t2,p[2]-n2[2]*t2];
                                                                        let v0b = [b2[0]-a2[0],b2[1]-a2[1],b2[2]-a2[2]]; let v1b = [c2[0]-a2[0],c2[1]-a2[1],c2[2]-a2[2]]; let v2b = [q2[0]-a2[0],q2[1]-a2[1],q2[2]-a2[2]];
                                                                        let d00b=v0b[0]*v0b[0]+v0b[1]*v0b[1]+v0b[2]*v0b[2]; let d01b=v0b[0]*v1b[0]+v0b[1]*v1b[1]+v0b[2]*v1b[2]; let d11b=v1b[0]*v1b[0]+v1b[1]*v1b[1]+v1b[2]*v1b[2]; let d20b=v2b[0]*v0b[0]+v2b[1]*v0b[1]+v2b[2]*v0b[2]; let d21b=v2b[0]*v1b[0]+v2b[1]*v1b[1]+v2b[2]*v1b[2];
                                                                        let denb=(d00b*d11b-d01b*d01b).max(1e-12); let wbb=(d11b*d20b-d01b*d21b)/denb; let wcc=(d00b*d21b-d01b*d20b)/denb; let waa=1.0-wbb-wcc;
                                                                        let mut kneg_c = 0usize; let mut vminc = waa;
                                                                        if wbb < vminc { vminc = wbb; kneg_c = 1; }
                                                                        if wcc < vminc { kneg_c = 2; }
                                                                        if waa < -1e-6 || wbb < -1e-6 || wcc < -1e-6 { let nf = neighbors[f_cpu as usize][kneg_c]; if nf!=u32::MAX { f_cpu = nf; } }
                                                                        let _ = writeln!(fcsv, "{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}", x, y, face_gpu, f_sim, f_cpu, kneg_s, kneg_c, wa, wb, wc, waa, wbb, wcc);
                                                                    }}
                                                                    println!("[parity] wrote parity_roi.csv");
                                                                }
                                                                // Also dump only mismatching pixels across the full frame (capped)
                                                                if let Ok(mut fmm) = std::fs::File::create("parity_mismatch.csv") {
                                                                    use std::io::Write;
                                                                    let mut count: usize = 0;
                                                                    let picker = GeoPicker::new();
                                                                    'outer: for y in 0..h { for x in 0..w {
                                                                        let idx = ((y * w + x) * 2) as usize; let face_gpu = dbg_buf[idx];
                                                                        let fx = x as f32 + 0.5; let fy = y as f32 + 0.5;
                                                                        let lon = (fx / w as f32) * std::f32::consts::TAU - std::f32::consts::PI;
                                                                        let lat = std::f32::consts::FRAC_PI_2 - (fy / h as f32) * std::f32::consts::PI;
                                                                        let cl = lat.cos(); let p = Vec3::new(cl*lon.cos(), lat.sin(), -cl*lon.sin()).norm();
                                                                        let cpu = picker.pick_from_unit(p, f);
                                                                        if face_gpu != cpu.face {
                                                                            let _ = writeln!(fmm, "{},{},{},{},{},{:.6},{:.6},{:.6}", x, y, face_gpu, cpu.face, cpu.kneg, cpu.w[0], cpu.w[1], cpu.w[2]);
                                                                            count += 1; if count >= 20000 { break 'outer; }
                                                                        }
                                                                    }}
                                                                    println!("[parity] wrote parity_mismatch.csv ({} rows)", 20000);
                                                                }
                                                            }
                                                        }
                                                    }
                                                    // Deferred CSV export (same pass) when requested
                                                    if ov.export_parity_csv_requested {
                                                        if let Some(buf) = rg.read_debug_face_tri(&gpu.device, &gpu.queue) {
                                                            // Export existing parity view for reference
                                                            let path = std::path::Path::new("parity_debug.csv");
                                                            if let Ok(mut file) = std::fs::File::create(path) {
                                                                use std::io::Write;
                                                                let _ = writeln!(file, "x,y,lon,lat,face_gpu,tri_gpu_global,face_cpu,tri_cpu_global,diff");
                                                                let (rw, rh) = ov.raster_size; let f = f_now; let picker = GeoPicker::new();
                                                                for y in 0..rh { for x in 0..rw {
                                                                    let idx=((y*rw+x)*2) as usize; let face_gpu=buf[idx]; let tri_gpu_global=buf[idx+1];
                                                                    let fx=x as f32+0.5; let fy=y as f32+0.5; let lon=(fx / rw as f32)*std::f32::consts::TAU - std::f32::consts::PI; let lat=std::f32::consts::FRAC_PI_2 - (fy / rh as f32)*std::f32::consts::PI; let cl=lat.cos(); let p=Vec3::new(cl*lon.cos(), lat.sin(), -cl*lon.sin()).norm();
                                                                    let cpu = picker.pick_from_unit(p, f);
                                                                    let tri_cpu_global = cpu.face * f * f + cpu.tri; // unify to global scheme like WGSL
                                                                    let diff=(face_gpu!=cpu.face)||(tri_gpu_global!=tri_cpu_global);
                                                                    let _=writeln!(file, "{},{},{:.6},{:.6},{},{},{},{},{}", x, y, lon, lat, face_gpu, tri_gpu_global, cpu.face, tri_cpu_global, if diff {1} else {0});
                                                                } }
                                                            }
                                                            println!("[parity] wrote parity_debug.csv ({}x{})", ov.raster_size.0, ov.raster_size.1);

                                                            // Export rollover probe CSV per checklist
                                                            let path_probe = std::path::Path::new("rollover_probe.csv");
                                                            if let Ok(mut file2) = std::fs::File::create(path_probe) {
                                                                use std::io::Write;
                                                                let _ = writeln!(file2, "x,y,f0,kneg,nf,f1,face_gpu,tri_gpu,wa0,wb0,wc0,wa1,wb1,wc1");
                                                                let (rw, rh) = ov.raster_size; let picker = GeoPicker::new();
                                                                // Access canonical faces from the shared picker to compute f0/kneg/nf
                                                                let faces = picker.faces();
                                                                for y in 0..rh { for x in 0..rw {
                                                                    let idx = ((y*rw + x) * 2) as usize;
                                                                    let face_gpu = buf[idx];
                                                                    let tri_gpu  = buf[idx + 1];
                                                                    // Pixel center mapping must match WGSL exactly
                                                                    let fx = x as f32 + 0.5; let fy = y as f32 + 0.5;
                                                                    let lon = (fx / rw as f32) * std::f32::consts::TAU - std::f32::consts::PI;
                                                                    let lat = std::f32::consts::FRAC_PI_2 - (fy / rh as f32) * std::f32::consts::PI;
                                                                    let cl = lat.cos(); let p = aule_geo::Vec3 { x: cl*lon.cos(), y: lat.sin(), z: -cl*lon.sin() };
                                                                    // f0: argmax before rollover
                                                                    let f0 = aule_geo::pick_face(p, faces);
                                                                    // Pre-rollover planar barycentrics on f0 (wa0,wb0,wc0)
                                                                    let w0 = aule_geo::barycentrics_plane(p, &faces[f0 as usize]);
                                                                    // kneg: index of most-negative component
                                                                    let mut kneg: u32 = 0; let mut minv = w0[0];
                                                                    if w0[1] < minv { kneg = 1; minv = w0[1]; }
                                                                    if w0[2] < minv { kneg = 2; /* minv = w0[2]; */ }
                                                                    // nf: neighbor across edge opposite kneg (table lookup)
                                                                    let nf = faces[f0 as usize].neighbor_opp[kneg as usize];
                                                                    // f1 per shader path: use GPU-written face as ground truth
                                                                    let f1 = face_gpu;
                                                                    // Post-rollover barycentrics on f1 (wa1,wb1,wc1)
                                                                    let w1 = aule_geo::barycentrics_plane(p, &faces[f1 as usize]);
                                                                    let _ = writeln!(file2, "{},{},{},{},{},{},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}", x, y, f0, kneg, nf, f1, face_gpu, tri_gpu, w0[0], w0[1], w0[2], w1[0], w1[1], w1[2]);
                                                                } }
                                                            }
                                                            println!("[probe] wrote rollover_probe.csv ({}x{})", ov.raster_size.0, ov.raster_size.1);
                                                        }
                                                        ov.export_parity_csv_requested = false;
                                                    }
                                                    if ov.raster_tex_id.is_none() { let tid = egui_renderer.register_native_texture(&gpu.device, &rg.out_view, wgpu::FilterMode::Linear); ov.raster_tex_id = Some(tid); }
                                                    ov.raster_dirty = false; ov.world_dirty = false; ov.color_dirty = false; ov.last_raster_at = std::time::Instant::now();
                                                    println!("[viewer] raster(gpu) W={} H={} | F={} | verts={} | face_tbl={} | dispatch={}x{}", rw, rh, f_now, world.depth_m.len(), 20 * ((f_now + 1) * (f_now + 2) / 2), (rw + 7) / 8, (rh + 7) / 8);
                                                }
                                            }
                                            if let Some(tid) = ov.raster_tex_id { let avail = ui.available_size(); ui.image(egui::load::SizedTexture::new(tid, avail)); }
                                            // Draw parity heat overlay as red dots if enabled
                                            if ov.show_parity_heat {
                                                if let Some(points) = &ov.parity_points {
                                                    let rect_px = rect;
                                                    // current raster size
                                                    let (rw, rh) = ov.raster_size;
                                                    let mut mesh = egui::epaint::Mesh::default();
                                                    for pxy in points.iter() {
                                                        let x = rect_px.left() + (pxy[0] as f32 + 0.5) / (rw as f32) * rect_px.width();
                                                        let y = rect_px.top() + (pxy[1] as f32 + 0.5) / (rh as f32) * rect_px.height();
                                                        overlay::OverlayState::mesh_add_dot(&mut mesh, egui::pos2(x, y), 1.2, egui::Color32::RED);
                                                    }
                                                    painter.add(egui::Shape::mesh(mesh));
                                                }
                                            }
                                            // Draw overlays on top (skip bathy points to avoid overdraw)
                                            let saved = ov.show_bathy; ov.show_bathy = false; overlay::draw_advanced_layers(ui, &painter, rect, &world, &world.grid, &mut ov); ov.show_bathy = saved;
                                        }
                                    } else {
                                        // CPU fallback raster
                                        if ov.raster_dirty || ov.raster_tex.is_none() {
                                            let (rw, rh) = ov.raster_size;
                                            let img = raster::render_map(&world, &ov, rw, rh);
                                            ov.raster_tex = Some(ctx.load_texture(
                                                "raster",
                                                img,
                                                egui::TextureOptions::LINEAR,
                                            ));
                                            ov.raster_dirty = false;
                                        }
                                        if let Some(tex) = &ov.raster_tex {
                                            painter.image(
                                                tex.id(),
                                                rect,
                                                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                                                egui::Color32::WHITE,
                                            );
                                        }
                                    }
                                    // Simple per-frame stepping using centralized pipeline (T-632)
                                    if ov.run_active && world.clock.t_myr < ov.run_target_myr {
                                        let n = ov.steps_per_frame.max(1);
                                        let mut steps_done = 0u32;
                                        for _ in 0..n {
                                            if world.clock.t_myr >= ov.run_target_myr { break; }
                                            let cfg = engine::pipeline::PipelineCfg {
                                                dt_myr: 1.0,
                                                steps_per_frame: 1,
                                                enable_flexure: !ov.disable_flexure,
                                                enable_erosion: !ov.disable_erosion,
                                                target_land_frac: ov.simple_target_land,
                                                freeze_eta: ov.freeze_eta,
                                                log_mass_budget: false,
                                                enable_subduction: !ov.disable_subduction,
                                                enable_rigid_motion: true,
                                                cadence_trf_every: ov.cadence_trf_every.max(1),
                                                cadence_sub_every: ov.cadence_sub_every.max(1),
                                                cadence_flx_every: ov.cadence_flx_every.max(1),
                                                cadence_sea_every: ov.cadence_sea_every.max(1),
                                                cadence_surf_every: ov.cadence_sea_every.max(1),
                                                substeps_transforms: 4,
                                                substeps_subduction: 4,
                                                use_gpu_flexure: false,
                                                gpu_flex_levels: ov.levels.max(1),
                                                gpu_flex_cycles: ov.flex_cycles.max(1),
                                                gpu_wj_omega: ov.wj_omega,
                                                subtract_mean_load: ov.subtract_mean_load,
                                                surf_k_stream: ov.surf_k_stream,
                                                surf_m_exp: ov.surf_m_exp,
                                                surf_n_exp: ov.surf_n_exp,
                                                surf_k_diff: ov.surf_k_diff,
                                                surf_k_tr: ov.surf_k_tr,
                                                surf_p_exp: ov.surf_p_exp,
                                                surf_q_exp: ov.surf_q_exp,
                                                surf_rho_sed: ov.surf_rho_sed,
                                                surf_min_slope: ov.surf_min_slope,
                                                surf_subcycles: ov.surf_subcycles.max(1),
                                                surf_couple_flexure: ov.surf_couple_flexure,
                                                sub_tau_conv_m_per_yr: ov.sub_tau_conv_m_per_yr,
                                                sub_trench_half_width_km: ov.sub_trench_half_width_km,
                                                sub_arc_offset_km: ov.sub_arc_offset_km,
                                                sub_arc_half_width_km: ov.sub_arc_half_width_km,
                                                sub_backarc_width_km: ov.sub_backarc_width_km,
                                                sub_trench_deepen_m: ov.sub_trench_deepen_m,
                                                sub_arc_uplift_m: ov.sub_arc_uplift_m,
                                                sub_backarc_uplift_m: ov.sub_backarc_uplift_m,
                                                sub_rollback_offset_m: ov.sub_rollback_offset_m,
                                                sub_rollback_rate_km_per_myr: ov.sub_rollback_rate_km_per_myr,
                                                sub_backarc_extension_mode: ov.sub_backarc_extension_mode,
                                                sub_backarc_extension_deepen_m: ov.sub_backarc_extension_deepen_m,
                                                sub_continent_c_min: ov.sub_continent_c_min,
                                                cadence_spawn_plate_every: 0,
                                                cadence_retire_plate_every: 0,
                                                cadence_force_balance_every: 8,
                                                fb_gain: 1.0e-12,
                                                fb_damp_per_myr: 0.2,
                                                fb_k_conv: 1.0,
                                                fb_k_div: 0.5,
                                                fb_k_trans: 0.1,
                                                fb_max_domega: 5.0e-9,
                                                fb_max_omega: 2.0e-7,
                                            };
                                            // Use world.depth_m order as backing for elevation upload later
                                            let mut elev_tmp: Vec<f32> = vec![0.0; world.depth_m.len()];
                                            let mut eta = world.sea.eta_m;
                                            engine::pipeline::step_full(
                                                &mut world,
                                                engine::pipeline::SurfaceFields { elev_m: &mut elev_tmp, eta_m: &mut eta },
                                                cfg,
                                            );
                                            world.sea.eta_m = eta;
                                            // Elevation stats and logging: z in elev_tmp is already (η − depth); do not add η again
                                            // Invalidate plate/boundary caches so overlays reflect updated plate_id and edges
                                            ov.plates_cache = None;
                                            ov.bounds_cache = None;
                                            let mut zmin = f32::INFINITY;
                                            let mut zmax = f32::NEG_INFINITY;
                                            let mut zsum = 0.0f64;
                                            let mut land_area = 0.0f64;
                                            let mut total_area = 0.0f64;
                                            for (&z0, &a) in elev_tmp.iter().zip(world.area_m2.iter()) {
                                                 if z0.is_finite() {
                                                     let z = z0; // already η − depth
                                                     zmin = zmin.min(z);
                                                     zmax = zmax.max(z);
                                                     zsum += z as f64;
                                                     total_area += a as f64;
                                                     if z > 0.0 { land_area += a as f64; }
                                                 }
                                             }
                                            let zmean = if total_area > 0.0 { zsum / total_area } else { 0.0 };
                                            println!("[draw] elev min/mean/max = {:.0}/{:.0}/{:.0} m | land={:.1}%", zmin, zmean, zmax, 100.0 * land_area / total_area.max(1.0));
                                            // Bookkeeping
                                            steps_done += 1;
                                            ov.color_dirty = true;
                                            ov.world_dirty = true;
                                            ctx.request_repaint();
                                            // Per-step land fraction (area-weighted) using elev = η − depth
                                            let area_total: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                                            let mut land_area: f64 = 0.0;
                                            for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) {
                                                if (world.sea.eta_m as f64 - d as f64) > 0.0 { land_area += a as f64; }
                                            }
                                            let land_frac = if area_total > 0.0 { land_area / area_total } else { 0.0 };
                                            println!(
                                                "[step/ui] Simple: land={:.1}% | eta={:+.0} m | steps/frame={}",
                                                land_frac * 100.0,
                                                world.sea.eta_m,
                                                n
                                            );
                                        }
                                        if steps_done > 0 && ov.debug_burst_steps > 0 { ov.debug_burst_steps = ov.debug_burst_steps.saturating_sub(steps_done); if ov.debug_burst_steps == 0 { println!("[debug] profiling burst complete."); } }
                                        if world.clock.t_myr >= ov.run_target_myr { ov.run_active = false; }
                                    }
                                } else {
                                    overlay::draw_advanced_layers(ui, &painter, rect, &world, &world.grid, &mut ov);
                                }
                            });
                            // Remove legacy central/top UI panels (moved to drawer)
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
                        // If globe path is enabled, render globe before egui
                        if ov.view_mode == overlay::ViewMode::Globe {
                            if let (Some(gr), Some(mesh)) = (pipes.globe.as_ref(), pipes.globe_mesh.as_ref()) {
                                // Update uniforms and draw
                                // Apply UI exaggeration to uniforms
                                let view_proj = if let Some(cam) = pipes.globe_cam.as_ref() { cam.view_proj() } else { glam::Mat4::IDENTITY.to_cols_array_2d() };
                                let radius = 1.0f32;
                                let exagger = ov.globe.exaggeration;
                                let dbg_flags = 0u32;
                                // Keep LUT in sync with overlay palette
                                gr.write_lut_from_overlay(&gpu.queue, &ov);
                                // If world changed, refresh vertex heights
                                if ov.world_dirty { gr.upload_heights(&gpu.queue, &world.depth_m); }
                                gr.update_uniforms(
                                    &gpu.queue,
                                    view_proj,
                                    radius,
                                    exagger,
                                    dbg_flags,
                                    ov.hypso_d_max,
                                    ov.hypso_h_max,
                                );
                                let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                                    label: Some("globe pass"),
                                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                                        view: &view,
                                        resolve_target: None,
                                        ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color { r: 0.02, g: 0.02, b: 0.04, a: 1.0 }), store: wgpu::StoreOp::Store },
                                    })],
                                    depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                                        view: &gpu.depth_view,
                                        depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                                        stencil_ops: None,
                                    }),
                                    occlusion_query_set: None,
                                    timestamp_writes: None,
                                });
                                gr.draw(&mut rpass, mesh);
                                if ov.globe.show_wireframe { gr.draw_lines(&mut rpass, mesh); }
                            }
                        }

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
                            let max_steps_frame = 8u32; // clamp to keep FPS ~60
                            while step_accum_s >= step_interval && steps_this_frame < max_steps_frame {
                                let _burst = ov.debug_burst_steps > 0;
                                let _enable_all = ov.debug_enable_all || _burst;
                                let _adv_every = ov.cadence_adv_every.max(1);
                                let _trf_every = ov.cadence_trf_every.max(1);
                                let _sub_every = ov.cadence_sub_every.max(1);
                                let _flx_every = ov.cadence_flx_every.max(1);
                                let _sea_every = ov.cadence_sea_every.max(1);
                                // Unified pipeline stepping path
                                let t0 = world.clock.t_myr;
                                let cfg = engine::pipeline::PipelineCfg {
                                    dt_myr,
                                    steps_per_frame: 1,
                                    enable_flexure: !ov.disable_flexure,
                                    enable_erosion: !ov.disable_erosion,
                                    target_land_frac: ov.simple_target_land,
                                    freeze_eta: ov.freeze_eta,
                                    log_mass_budget: false,
                                    enable_subduction: !ov.disable_subduction,
                                    enable_rigid_motion: true,
                                    cadence_trf_every: ov.cadence_trf_every.max(1),
                                    cadence_sub_every: ov.cadence_sub_every.max(1),
                                    cadence_flx_every: ov.cadence_flx_every.max(1),
                                    cadence_sea_every: ov.cadence_sea_every.max(1),
                                    cadence_surf_every: ov.cadence_sea_every.max(1),
                                    substeps_transforms: 4,
                                    substeps_subduction: 4,
                                    use_gpu_flexure: false,
                                    gpu_flex_levels: ov.levels.max(1),
                                    gpu_flex_cycles: ov.flex_cycles.max(1),
                                    gpu_wj_omega: ov.wj_omega,
                                    subtract_mean_load: ov.subtract_mean_load,
                                    surf_k_stream: ov.surf_k_stream,
                                    surf_m_exp: ov.surf_m_exp,
                                    surf_n_exp: ov.surf_n_exp,
                                    surf_k_diff: ov.surf_k_diff,
                                    surf_k_tr: ov.surf_k_tr,
                                    surf_p_exp: ov.surf_p_exp,
                                    surf_q_exp: ov.surf_q_exp,
                                    surf_rho_sed: ov.surf_rho_sed,
                                    surf_min_slope: ov.surf_min_slope,
                                    surf_subcycles: ov.surf_subcycles.max(1),
                                    surf_couple_flexure: ov.surf_couple_flexure,
                                    sub_tau_conv_m_per_yr: ov.sub_tau_conv_m_per_yr,
                                    sub_trench_half_width_km: ov.sub_trench_half_width_km,
                                    sub_arc_offset_km: ov.sub_arc_offset_km,
                                    sub_arc_half_width_km: ov.sub_arc_half_width_km,
                                    sub_backarc_width_km: ov.sub_backarc_width_km,
                                    sub_trench_deepen_m: ov.sub_trench_deepen_m,
                                    sub_arc_uplift_m: ov.sub_arc_uplift_m,
                                    sub_backarc_uplift_m: ov.sub_backarc_uplift_m,
                                    sub_rollback_offset_m: ov.sub_rollback_offset_m,
                                    sub_rollback_rate_km_per_myr: ov.sub_rollback_rate_km_per_myr,
                                    sub_backarc_extension_mode: ov.sub_backarc_extension_mode,
                                    sub_backarc_extension_deepen_m: ov.sub_backarc_extension_deepen_m,
                                    sub_continent_c_min: ov.sub_continent_c_min,
                                    cadence_spawn_plate_every: 0,
                                    cadence_retire_plate_every: 0,
                                    cadence_force_balance_every: 8,
                                    fb_gain: 1.0e-12,
                                    fb_damp_per_myr: 0.2,
                                    fb_k_conv: 1.0,
                                    fb_k_div: 0.5,
                                    fb_k_trans: 0.1,
                                    fb_max_domega: 5.0e-9,
                                    fb_max_omega: 2.0e-7,
                                };
                                let mut elev_tmp: Vec<f32> = vec![0.0; world.depth_m.len()];
                                let mut eta = world.sea.eta_m;
                                engine::pipeline::step_full(
                                    &mut world,
                                    engine::pipeline::SurfaceFields { elev_m: &mut elev_tmp, eta_m: &mut eta },
                                    cfg,
                                );
                                world.sea.eta_m = eta;
                                // Step log (one line per step) using elev = η − depth
                                let area_total: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                                let mut land_area: f64 = 0.0;
                                for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) {
                                    if (world.sea.eta_m as f64 - d as f64) > 0.0 { land_area += a as f64; }
                                }
                                let land_frac = if area_total > 0.0 { land_area / area_total } else { 0.0 };
                                println!(
                                    "[step] t={:.1}→{:.1} Myr | land={:.1}% | residual flex={}",
                                    t0, world.clock.t_myr, land_frac * 100.0,
                                    if !ov.disable_flexure { format!("{:.3e}", world.last_flex_residual) } else { "–".to_string() }
                                );
                                // Edge-triggered snapshots
                                let f = snapshot_every_myr.max(0.0);
                                if f > 0.0 {
                                    if next_snapshot_t.is_infinite() || next_snapshot_t.is_nan() {
                                        let k = (world.clock.t_myr / (f as f64)).ceil();
                                        next_snapshot_t = (k * (f as f64)).max(world.clock.t_myr);
                                    }
                                    if world.clock.t_myr + 1e-9 >= next_snapshot_t {
                                        let name = format!("depth_t{:08.1}Myr.csv", world.clock.t_myr);
                                        let path = std::path::Path::new(&name);
                                        let _ = engine::snapshots::write_csv_depth(path, world.clock.t_myr, &world.depth_m);
                                        println!("[snapshot] wrote {} (N={})", name, world.depth_m.len());
                                        next_snapshot_t += f as f64;
                                    }
                                } else {
                                    next_snapshot_t = f64::INFINITY;
                                }
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

// Helper to apply preset values to subduction controls
fn apply_sub_preset(ov: &mut overlay::OverlayState) {
    match ov.sub_preset {
        1 => {
            // Strong rollback
            ov.sub_rollback_offset_m = 80_000.0;
            ov.sub_rollback_rate_km_per_myr = 20.0;
            ov.sub_trench_deepen_m = 2200.0;
            ov.sub_arc_uplift_m = -400.0;
            ov.sub_backarc_uplift_m = -150.0;
            ov.sub_backarc_extension_mode = false;
        }
        2 => {
            // Back-arc extension
            ov.sub_backarc_extension_mode = true;
            ov.sub_backarc_extension_deepen_m = 600.0;
            ov.sub_trench_deepen_m = 2000.0;
            ov.sub_arc_uplift_m = -250.0;
            ov.sub_backarc_uplift_m = 0.0;
            ov.sub_rollback_offset_m = 40_000.0;
            ov.sub_rollback_rate_km_per_myr = 10.0;
        }
        3 => {
            // Weak arcs
            ov.sub_trench_deepen_m = 1500.0;
            ov.sub_arc_uplift_m = -120.0;
            ov.sub_backarc_uplift_m = -60.0;
            ov.sub_backarc_extension_mode = false;
            ov.sub_rollback_offset_m = 0.0;
            ov.sub_rollback_rate_km_per_myr = 0.0;
        }
        _ => {
            // Reference
            ov.sub_tau_conv_m_per_yr = 0.005;
            ov.sub_trench_half_width_km = 40.0;
            ov.sub_arc_offset_km = 140.0;
            ov.sub_arc_half_width_km = 25.0;
            ov.sub_backarc_width_km = 120.0;
            ov.sub_trench_deepen_m = 1800.0;
            ov.sub_arc_uplift_m = -300.0;
            ov.sub_backarc_uplift_m = -120.0;
            ov.sub_rollback_offset_m = 0.0;
            ov.sub_rollback_rate_km_per_myr = 0.0;
            ov.sub_backarc_extension_mode = false;
            ov.sub_backarc_extension_deepen_m = 400.0;
            ov.sub_continent_c_min = 0.6;
        }
    }
}
