//! Aulë viewer binary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

mod colormap;
mod overlay;
mod plot;
mod plot_age_depth;
mod plot_flexure;

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
            let _ = engine::world::step_once(world, sp);
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
                ctx.request_repaint();
            }
        }
        ui.add(egui::DragValue::new(&mut ov.simple_t_end_myr).speed(10.0).suffix(" Myr"));

        if ui.button("Generate world").clicked() {
            // a) Reset world (keep F)
            let f = ov.simple_f.max(1);
            let (preset_n_continents, preset_radius_km, preset_falloff_km, preset_plates): (
                u32,
                f64,
                f64,
                u32,
            ) = match ov.simple_preset {
                1 => (4, 2600.0, 700.0, 10),
                2 => (2, 1800.0, 500.0, 6),
                _ => (3, 2200.0, 600.0, 8),
            };
            *world = engine::world::World::new(f, preset_plates, ov.simple_seed);
            // b) Apply preset continents to target land fraction
            let cp = engine::continent::ContinentParams {
                seed: ov.simple_seed,
                n_continents: preset_n_continents,
                mean_radius_km: preset_radius_km,
                falloff_km: preset_falloff_km,
                plateau_uplift_m: 1.0,
                target_land_fraction: None,
            };
            let cf = engine::continent::build_continents(&world.grid, cp);
            // Baseline from age only first
            for (i, d) in world.depth_m.iter_mut().enumerate() {
                let mut v =
                    engine::age::depth_from_age(world.age_myr[i] as f64, 2600.0, 350.0, 0.0) as f32;
                if !v.is_finite() {
                    v = 6000.0;
                }
                *d = v.clamp(0.0, 6000.0);
            }
            // Solve amplitude for target land fraction and apply
            let amp_m = engine::continent::solve_amplitude_for_target_land_fraction(
                &world.depth_m,
                &cf.uplift_template_m,
                &world.area_m2,
                ov.simple_target_land as f64,
                1e-3,
                64,
            );
            let (_mask_land, _land_frac) = engine::continent::apply_continents(
                &mut world.depth_m,
                &cf.uplift_template_m,
                amp_m,
                &world.area_m2,
            );
            // Diagnostics after uplift
            let mut dmin = f32::INFINITY;
            let mut dmax = f32::NEG_INFINITY;
            let mut dsum = 0.0f64;
            for &d in &world.depth_m {
                if d.is_finite() {
                    if d < dmin {
                        dmin = d;
                    }
                    if d > dmax {
                        dmax = d;
                    }
                    dsum += d as f64;
                }
            }
            let dmean = dsum / (world.depth_m.len().max(1) as f64);
            println!("[cont] amp={:.0} m  land_target={:.1}%  after_uplift depth min/mean/max={:.0}/{:.0}/{:.0} m",
                amp_m, ov.simple_target_land * 100.0, dmin, dmean, dmax);
            // c) Apply palette before running
            apply_simple_palette(ov, world, ctx);
            // d) Solve sea level offset to hit target ocean fraction = 1 - target_land
            // Capture reference volume/area if missing
            if world.sea_level_ref.is_none() {
                world.sea_level_ref = Some(engine::isostasy::compute_ref(&world.depth_m, &world.area_m2));
            }
            if let Some(ref_ref) = world.sea_level_ref {
                let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                let frac_ref = if total_area > 0.0 { ref_ref.ocean_area_m2 / total_area } else { 0.0 };
                let target_ocean = 1.0 - (ov.simple_target_land as f64);
                let target_volume = if frac_ref > 0.0 { ref_ref.volume_m3 * (target_ocean / frac_ref) } else { ref_ref.volume_m3 };
                let off = engine::isostasy::solve_offset_for_volume(&world.depth_m, &world.area_m2, target_volume, 1e6, 64);
                for d in &mut world.depth_m { *d = (*d as f64 + off) as f32; }
                let ocean_cells = world.depth_m.iter().filter(|d| **d > 0.0).count();
                let ocean_frac = (ocean_cells as f64) / (world.depth_m.len().max(1) as f64);
                println!("[sea]  policy=target  ocean={:.1}%  offset={:.0} m", ocean_frac * 100.0, off);
            }
            // e) Run to t_end (safe dt) with realtime redraw
            let sp = engine::world::StepParams {
                dt_myr: 1.0,
                do_flexure: false,
                do_isostasy: true,
                do_transforms: false,
                do_subduction: false,
                do_continents: false,
                do_ridge_birth: true,
                auto_rebaseline_after_continents: false,
                do_rigid_motion: true,
                do_orogeny: false,
                do_accretion: false,
                do_rifting: false,
                do_surface: false,
                surface_params: engine::surface::SurfaceParams {
                    k_stream: ov.surf_k_stream,
                    m_exp: ov.surf_m_exp,
                    n_exp: ov.surf_n_exp,
                    k_diff: ov.surf_k_diff,
                    k_tr: ov.surf_k_tr,
                    p_exp: ov.surf_p_exp,
                    q_exp: ov.surf_q_exp,
                    rho_sed: ov.surf_rho_sed,
                    min_slope: ov.surf_min_slope,
                    subcycles: ov.surf_subcycles.max(1),
                    couple_flexure: ov.surf_couple_flexure,
                },
            };
            run_to_t_realtime(ctx, world, ov, &sp, ov.simple_t_end_myr, 4);
            // f) Start realtime run; per-frame updater handles drawing
            ov.run_active = true;
            ov.run_target_myr = ov.simple_t_end_myr;
            ctx.request_repaint();
        }
    });

    egui::CollapsingHeader::new("Continents & Seeds").default_open(true).show(ui, |ui| {
        ui.horizontal(|ui| {
            ui.label("Preset:");
            // dropdown for ov.simple_preset (WorldPreset)
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

    egui::CollapsingHeader::new("Surface").default_open(ov.adv_open_surface).show(ui, |ui| {
        ui.checkbox(&mut ov.surface_enable, "Enable surface processes");
        ui.checkbox(&mut ov.surf_couple_flexure, "Couple to flexure");
        ui.add(
            egui::Slider::new(&mut ov.surf_k_stream, 1.0e-6..=1.0e-4)
                .logarithmic(true)
                .text("K_stream"),
        );
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut ov.surf_m_exp, 0.3..=0.6).text("m"));
            ui.add(egui::Slider::new(&mut ov.surf_n_exp, 0.8..=1.5).text("n"));
        });
        ui.add(
            egui::Slider::new(&mut ov.surf_k_diff, 0.05..=0.5)
                .logarithmic(true)
                .text("κ_diff (m²/yr)"),
        );
        ui.add(
            egui::Slider::new(&mut ov.surf_k_tr, 0.01..=0.3).logarithmic(true).text("K_transport"),
        );
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut ov.surf_p_exp, 1.0..=2.0).text("p"));
            ui.add(egui::Slider::new(&mut ov.surf_q_exp, 0.5..=1.5).text("q"));
        });
        ui.add(egui::Slider::new(&mut ov.surf_rho_sed, 1000.0..=2500.0).text("ρ_sed (kg/m³)"));
        ui.add(
            egui::Slider::new(&mut ov.surf_min_slope, 1.0e-5..=1.0e-3)
                .logarithmic(true)
                .text("min slope"),
        );
        let mut sc = ov.surf_subcycles;
        let changed_sc = ui.add(egui::Slider::new(&mut sc, 1..=10).text("Subcycles")).changed();
        if changed_sc {
            ov.surf_subcycles = sc;
        }
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
                                    ui.selectable_value(&mut ov.mode_simple, true, "Simple");
                                    ui.selectable_value(&mut ov.mode_simple, false, "Advanced");
                                    if ui.button(if ov.drawer_open { "⟨⟩" } else { "☰" }).clicked() {
                                        ov.drawer_open = !ov.drawer_open;
                                    }
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
                                        });
                                });
                            }
                            // Central canvas (draw only, no controls)
                            egui::CentralPanel::default().show(ctx, |ui| {
                                let rect = ui.max_rect();
                                let painter = ui.painter_at(rect);
                                if ov.mode_simple {
                                    overlay::draw_color_layer(ui, &painter, rect, &world, &world.grid, &mut ov);
                                    // Simple per-frame stepping
                                    if ov.run_active && world.clock.t_myr < ov.run_target_myr {
                                        let n = ov.steps_per_frame.max(1);
                                        let mut steps_done = 0u32;
                                        for _ in 0..n {
                                            if world.clock.t_myr >= ov.run_target_myr { break; }
                                            let sp = engine::world::StepParams {
                                                dt_myr: 1.0,
                                                do_flexure: false,
                                                do_isostasy: true,
                                                do_transforms: false,
                                                do_subduction: false,
                                                do_continents: false,
                                                do_ridge_birth: true,
                                                auto_rebaseline_after_continents: false,
                                                do_rigid_motion: true,
                                                do_orogeny: false,
                                                do_accretion: false,
                                                do_rifting: false,
                                                do_surface: false,
                                                surface_params: engine::surface::SurfaceParams {
                                                    k_stream: ov.surf_k_stream,
                                                    m_exp: ov.surf_m_exp,
                                                    n_exp: ov.surf_n_exp,
                                                    k_diff: ov.surf_k_diff,
                                                    k_tr: ov.surf_k_tr,
                                                    p_exp: ov.surf_p_exp,
                                                    q_exp: ov.surf_q_exp,
                                                    rho_sed: ov.surf_rho_sed,
                                                    min_slope: ov.surf_min_slope,
                                                    subcycles: ov.surf_subcycles.max(1),
                                                    couple_flexure: ov.surf_couple_flexure,
                                                },
                                            };
                                            let t0 = world.clock.t_myr;
                                            let _ = engine::world::step_once(&mut world, &sp);
                                            steps_done += 1;
                                            ov.color_dirty = true;
                                            ov.world_dirty = true;
                                            ctx.request_repaint();
                                            println!("[step/ui] t={:.1}→{:.1} Myr (+{} steps)", t0, world.clock.t_myr, steps_done);
                                        }
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
                                let sp = engine::world::StepParams {
                                    dt_myr: dt_myr as f64,
                                    do_flexure: ov.enable_flexure,
                                    do_isostasy: ov.apply_sea_level,
                                    do_transforms: ov.show_transforms,
                                    do_subduction: ov.show_subduction,
                                    do_continents: ov.continents_apply,
                                    do_ridge_birth: true,
                                    auto_rebaseline_after_continents: ov.auto_rebaseline_l,
                                    do_rigid_motion: ov.kin_enable,
                                    do_orogeny: false,
                                    do_accretion: false,
                                    do_rifting: false,
                                    do_surface: ov.surface_enable,
                                    surface_params: engine::surface::SurfaceParams {
                                        k_stream: ov.surf_k_stream,
                                        m_exp: ov.surf_m_exp,
                                        n_exp: ov.surf_n_exp,
                                        k_diff: ov.surf_k_diff,
                                        k_tr: ov.surf_k_tr,
                                        p_exp: ov.surf_p_exp,
                                        q_exp: ov.surf_q_exp,
                                        rho_sed: ov.surf_rho_sed,
                                        min_slope: ov.surf_min_slope,
                                        subcycles: ov.surf_subcycles.max(1),
                                        couple_flexure: ov.surf_couple_flexure,
                                    },
                                };
                                let stats = engine::world::step_once(&mut world, &sp);
                                // Step log (one line per step)
                                let area_total: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
                                let land_area: f64 = world
                                    .depth_m
                                    .iter()
                                    .zip(world.area_m2.iter())
                                    .filter(|(&d, _)| d <= 0.0)
                                    .map(|(_, &a)| a as f64)
                                    .sum();
                                let land_frac = if area_total > 0.0 { land_area / area_total } else { 0.0 };
                                println!(
                                    "[step] t={:.1} Myr dt={:.1} | div={} conv={} trans={} | land={:.1}% C̄={:.1}% | residual flex={}",
                                    stats.t_myr, stats.dt_myr, stats.div_count, stats.conv_count, stats.trans_count, land_frac * 100.0, stats.c_bar * 100.0,
                                    if sp.do_flexure { format!("{:.3e}", world.last_flex_residual) } else { "–".to_string() }
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
