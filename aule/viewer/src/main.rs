//! Aulë viewer binary.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::dbg_macro, clippy::large_enum_variant)]

mod colormap;
mod globe;
mod gpu_buffer_manager;
mod overlay;
mod plot;
mod plot_age_depth;
mod plot_flexure;
mod raster;
mod raster_gpu;
use viewer::cpu_picker::{GeoPicker, Vec3};
use viewer::pixel_map::{pixel_to_lon_lat, sph_to_unit as sph_to_unit_px};
use gpu_buffer_manager::GpuBufferManager;


use egui_wgpu::Renderer as EguiRenderer;
use egui_wgpu::ScreenDescriptor;
use egui_winit::State as EguiWinitState;
use winit::{
    dpi::PhysicalSize,
    event::{Event, WindowEvent},
    event_loop::EventLoop,
    window::{Window, WindowBuilder},
};

// --- Hypsometry snapshot guards (poison detection & slew) ---
// Physical constants - use centralized values from engine
fn get_phys_consts() -> engine::PhysConsts {
    engine::PhysConsts::default()
}

// Removed unused imports

// --- Simulation thread messages and snapshot types ---
#[derive(Clone, Debug)]
struct WorldSnapshot {
    depth_m: Vec<f32>,
    c: Vec<f32>,
    th_c_m: Vec<f32>,
    sea_eta_m: f32,
    plate_id: Vec<u16>,
    pole_axis: Vec<[f32; 3]>,
    omega_rad_yr: Vec<f32>,
    v_en: Vec<[f32; 2]>,
    age_myr: Vec<f32>,
}

#[derive(Clone, Debug)]
struct ProcessFlags {
    pub enable_rigid_motion: bool,
    pub enable_subduction: bool,
    pub enable_transforms: bool,
    pub enable_flexure: bool,
    pub enable_surface_processes: bool,
    pub enable_isostasy: bool,
    pub enable_continental_buoyancy: bool,
    pub enable_orogeny: bool,
    pub enable_accretion: bool,
    pub enable_rifting: bool,
    pub enable_ridge_birth: bool,
    
    // Flexure backend configuration
    pub flexure_backend_cpu: bool,
    pub flexure_gpu_levels: u32,
    pub flexure_gpu_cycles: u32,
    
    // Unified cadence configuration
    pub cadence_config: engine::cadence_manager::CadenceConfig,
}

#[derive(Clone, Debug)]
enum SimCommand {
    Step(engine::config::PipelineCfg, ProcessFlags),
    SyncWorld(WorldSnapshot),
    Stop,
}

// T-505 drawer render — Simple/Advanced
#[allow(dead_code)]


/// Unified progressive UI that replaces Simple/Advanced modes
/// Features progressive disclosure with smart defaults and tooltips
fn render_unified_progressive_panels(
    ui: &mut egui::Ui,
    ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
    tx_cmd: &std::sync::mpsc::Sender<SimCommand>,
) {
    // 📋 Simulation Basics (Always expanded - essential controls)
    render_simulation_basics_unified(ui, ctx, world, ov, tx_cmd);
    
    // 🌍 World Generation (Expanded by default - common task)  
    render_world_generation_unified(ui, ctx, world, ov, tx_cmd);
    
    // 🎨 Visual Overlays (Collapsed by default - organized by complexity)
    render_visual_overlays_unified(ui, ctx, world, ov);
    
    // ⚙️ Physics Processes (Collapsed by default - advanced tuning)
    render_physics_processes_unified(ui, ctx, world, ov);
    
    // 🐛 Debug & Diagnostics (Collapsed by default - expert features)
    render_debug_diagnostics_unified(ui, ctx, world, ov);
}

/// 📋 Simulation Basics - Essential controls always visible
fn render_simulation_basics_unified(
    ui: &mut egui::Ui,
    ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
    _tx_cmd: &std::sync::mpsc::Sender<SimCommand>,
) {
    egui::CollapsingHeader::new("📋 Simulation Basics")
        .default_open(true)
        .show(ui, |ui| {
            // Resolution selector with availability indicators
            let tiling_ok = cfg!(feature = "tiling") || ov.high_f_available;
            let current_f = world.grid.frequency;
            ui.horizontal(|ui| {
                ui.label("Resolution (F):")
                    .on_hover_text("Grid resolution - higher values give more detail but slower simulation");
                ui.selectable_value(&mut ov.simple_f, 64, "64")
                    .on_hover_text("Fast - Good for learning and quick experiments");
                ui.selectable_value(&mut ov.simple_f, 128, "128")
                    .on_hover_text("Balanced - Recommended for most simulations");
                ui.selectable_value(&mut ov.simple_f, 256, "256")
                    .on_hover_text("High detail - Better for detailed analysis");
                ui.add_enabled_ui(tiling_ok, |ui| {
                    ui.selectable_value(&mut ov.simple_f, 512, "512")
                        .on_hover_text("Maximum detail - Requires high-end hardware");
                });
                if !tiling_ok {
                    ui.label(egui::RichText::new("(512 requires T-455)").small().color(egui::Color32::GRAY));
                }
            });
            
            // Rebuild world on resolution change
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
            
            ui.separator();
            
            // Time controls with smart defaults
            ui.horizontal(|ui| {
                ui.label("Time step:")
                    .on_hover_text("Simulation time step - smaller values are more accurate but slower");
                // Green highlight for recommended range (0.5-2.0 Myr)
                let mut dt = ov.sim_dt_myr;
                let slider = egui::Slider::new(&mut dt, 0.1..=5.0)
                    .suffix(" Myr")
                    .show_value(true);
                let response = ui.add(slider);
                
                // Highlight recommended range in green
                if dt >= 0.5 && dt <= 2.0 {
                    response.on_hover_text("✅ Recommended range for most simulations");
                } else if dt < 0.5 {
                    response.on_hover_text("⚠️ Very small - will be slow but very accurate");
                } else {
                    response.on_hover_text("⚠️ Large step - faster but may be less stable");
                }
                ov.sim_dt_myr = dt;
            });
            
            ui.horizontal(|ui| {
                ui.label("Target time:")
                    .on_hover_text("How long to run the simulation");
                ui.add(egui::DragValue::new(&mut ov.simple_t_end_myr)
                    .speed(10.0)
                    .suffix(" Myr")
                    .clamp_range(0.0..=2000.0));
            });
            
            ui.separator();
            
            // Simulation status and controls
            ui.horizontal(|ui| {
                ui.label(format!("Current: {:.1} Myr (Step {})", world.clock.t_myr, world.clock.step_idx));
                
                // Play/Pause button
                if ov.run_active {
                    if ui.button("⏸️ Pause")
                        .on_hover_text("Pause the simulation")
                        .clicked() {
                        ov.run_active = false;
                        ov.stepper.playing = false;
                    }
                    ui.label("🔄 Running");
                } else {
                    if ui.button("▶️ Play")
                        .on_hover_text("Start/resume the simulation")
                        .clicked() {
                        ov.run_active = true;
                        ov.stepper.playing = true;
                        ov.stepper.t_target_myr = ov.simple_t_end_myr as f32;
                    }
                    ui.label("⏸️ Paused");
                }
            });
        });
}

/// 🌍 World Generation - Progressive complexity from presets to detailed control
fn render_world_generation_unified(
    ui: &mut egui::Ui,
    ctx: &egui::Context,
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
    tx_cmd: &std::sync::mpsc::Sender<SimCommand>,
) {
    egui::CollapsingHeader::new("🌍 World Generation")
        .default_open(true)
        .show(ui, |ui| {
            // Quick preset buttons for beginners
            ui.label("Quick Presets:");
            ui.horizontal(|ui| {
                if ui.button("🏝️ Archipelago")
                    .on_hover_text("Many small islands - great for learning plate tectonics")
                    .clicked() {
                    ov.simple_preset = 1;
                    generate_world_with_preset(world, ov, ctx, tx_cmd);
                }
                if ui.button("🌍 Pangaea")
                    .on_hover_text("Supercontinent - shows continental breakup")
                    .clicked() {
                    ov.simple_preset = 3;
                    generate_world_with_preset(world, ov, ctx, tx_cmd);
                }
                if ui.button("🎲 Random")
                    .on_hover_text("Balanced random world - good default")
                    .clicked() {
                    ov.simple_preset = 0;
                    generate_world_with_preset(world, ov, ctx, tx_cmd);
                }
            });
            
            ui.separator();
            
            // Basic parameters with smart defaults
            ui.horizontal(|ui| {
                ui.label("Plates:")
                    .on_hover_text("Number of tectonic plates - more plates = more complex tectonics");
                let mut plates = match ov.simple_preset {
                    1 => 10,  // Archipelago
                    2 => 6,   // Continental
                    3 => 8,   // Supercontinent  
                    _ => 8,   // Default
                };
                let slider = egui::Slider::new(&mut plates, 4..=16).show_value(true);
                let response = ui.add(slider);
                
                // Highlight recommended range
                if plates >= 6 && plates <= 10 {
                    response.on_hover_text("✅ Good balance of complexity and performance");
                } else if plates < 6 {
                    response.on_hover_text("⚠️ Few plates - simpler but less realistic");
                } else {
                    response.on_hover_text("⚠️ Many plates - complex but slower simulation");
                }
                
                // Update preset if changed
                if plates != match ov.simple_preset { 1 => 10, 2 => 6, 3 => 8, _ => 8 } {
                    ov.simple_preset = 0; // Custom
                }
            });
            
            ui.horizontal(|ui| {
                ui.label("Land fraction:")
                    .on_hover_text("Percentage of surface that's land vs ocean");
                let slider = egui::Slider::new(&mut ov.simple_target_land, 0.1..=0.6)
                    .show_value(true)
                    .custom_formatter(|n, _| format!("{:.0}%", n * 100.0));
                let response = ui.add(slider);
                
                // Highlight Earth-like range  
                if ov.simple_target_land >= 0.25 && ov.simple_target_land <= 0.35 {
                    response.on_hover_text("✅ Earth-like land distribution");
                } else if ov.simple_target_land < 0.25 {
                    response.on_hover_text("🌊 Ocean world - mostly water");
                } else {
                    response.on_hover_text("🏔️ Continental world - mostly land");
                }
            });
            
            // Advanced options (collapsible)
            ui.collapsing("🔬 Advanced Options", |ui| {
                ui.horizontal(|ui| {
                    ui.label("Random seed:")
                        .on_hover_text("Controls random world generation - same seed = same world");
                    ui.add(egui::DragValue::new(&mut ov.simple_seed).speed(1000));
                    if ui.button("🎲").on_hover_text("Generate new random seed").clicked() {
                        ov.simple_seed = std::time::SystemTime::now()
                            .duration_since(std::time::UNIX_EPOCH)
                            .unwrap_or_default()
                            .as_secs();
                    }
                });
                
                ui.label("Preset details:");
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut ov.simple_preset, 0, "Balanced")
                        .on_hover_text("Good mix of land and ocean");
                    ui.selectable_value(&mut ov.simple_preset, 1, "Archipelago")
                        .on_hover_text("Many small continents");
                    ui.selectable_value(&mut ov.simple_preset, 2, "Continental")
                        .on_hover_text("Few large continents");
                    ui.selectable_value(&mut ov.simple_preset, 3, "Supercontinent")
                        .on_hover_text("One massive landmass");
                });
            });
            
            ui.separator();
            
            // Generate button
            if ui.button("🌍 Generate New World")
                .on_hover_text("Create a new world with current settings")
                .clicked() {
                generate_world_with_preset(world, ov, ctx, tx_cmd);
            }
        });
}

/// Helper function to generate world with current preset
fn generate_world_with_preset(
    world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
    ctx: &egui::Context,
    tx_cmd: &std::sync::mpsc::Sender<SimCommand>,
) {
    let plates = match ov.simple_preset {
        1 => 10,  // Archipelago
        2 => 6,   // Continental
        3 => 8,   // Supercontinent  
        _ => 8,   // Default
    };
    
    let continents_n = match ov.simple_preset {
        1 => 4,   // Archipelago - many small continents
        2 => 2,   // Continental - few large continents
        3 => 0,   // Supercontinent - triggers specialized generator
        _ => 3,   // Default - balanced
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
    
    // Classify boundaries and initialize
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
    
    // Set baseline oceanic bathymetry from age
    let n_cells = world.grid.cells;
    for i in 0..n_cells {
        let mut d = engine::age::depth_from_age_plate(
            world.age_myr[i] as f64, 
            world.plates.plate_id[i] as f64,
            0.0, 0.0, 0.0  // Additional required parameters
        );
        if d < 10.0 {
            d = 10.0; // Minimum depth
        }
        world.depth_m[i] = d as f32;
    }
    
    // Generate realistic continents with multifractal noise
    if continents_n > 0 {
        generate_realistic_continents(world, ov, continents_n);
        
        // FALLBACK: If no continental material was created, force create some
        let c_total_after: f32 = world.c.iter().sum();
        if c_total_after < 0.1 {
            println!("[world_gen] WARNING: No continental material created, applying fallback generation");
            create_simple_fallback_continents(world, ov, continents_n);
        }
    }
    
    // CRITICAL: Apply continental uplift to create visible land elevation
    println!("[world_gen] BEFORE uplift: depth_m range [{:.1}, {:.1}]", 
        world.depth_m.iter().fold(f32::INFINITY, |a, &b| a.min(b)),
        world.depth_m.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b)));
    
    let c_total: f32 = world.c.iter().sum();
    let th_total: f32 = world.th_c_m.iter().sum();
    println!("[world_gen] Continental data: C_sum={:.1}, th_c_sum={:.1} m", c_total, th_total);
    
    engine::continent::apply_uplift_from_c_thc(&mut world.depth_m, &world.c, &world.th_c_m);
    
    println!("[world_gen] AFTER uplift: depth_m range [{:.1}, {:.1}]", 
        world.depth_m.iter().fold(f32::INFINITY, |a, &b| a.min(b)),
        world.depth_m.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b)));
    
    // Add realistic topographic variation using multifractal noise
    add_realistic_topography(world, ov);
    
    // Mark elevation-dependent systems dirty
    ov.bathy_cache = None;
    ov.plates_cache = None;
    ov.vel_cache = None;
    ov.bounds_cache = None;
    ov.color_dirty = true;
    ov.world_dirty = true;
    ov.raster_dirty = true;
    // Re-enable GPU raster for smooth rendering after generation
    ov.use_gpu_raster = true;
    ov.raster_tex = None;
    ov.raster_tex_id = None;
    
    // GPU raster will provide smooth elevation colors; disable redundant CPU overlay
    // (The unified logic in main loop will set show_bathy = !use_gpu_raster)
    
    // CRITICAL: Sync the generated world to the simulation thread
    let ws = WorldSnapshot {
        depth_m: world.depth_m.clone(),
        c: world.c.clone(),
        th_c_m: world.th_c_m.clone(),
        sea_eta_m: world.sea.eta_m,
        plate_id: world.plates.plate_id.clone(),
        pole_axis: world.plates.pole_axis.clone(),
        omega_rad_yr: world.plates.omega_rad_yr.clone(),
        v_en: world.v_en.clone(),
        age_myr: world.age_myr.clone(),
    };
    let _ = tx_cmd.send(SimCommand::SyncWorld(ws));
    
    // Force repaint immediately
    ctx.request_repaint();
}

/// Generate realistic continents with multifractal noise and plate-based seeding
fn generate_realistic_continents(
    world: &mut engine::world::World,
    ov: &overlay::OverlayState,
    continents_n: u32,
) {
    let _n_cells = world.grid.cells;
    let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
    let target_land_area_m2 = ov.simple_target_land as f64 * total_area;
    
    // Clear existing continental data
    world.c.fill(0.0);
    world.th_c_m.fill(0.0);
    
    // Classify plates as continental or oceanic based on their properties
    let continental_plates = classify_continental_plates(world, ov);
    
    println!("[world_gen] Generating {} continents on {} continental plates", continents_n, continental_plates.len());
    
    // Generate continents preferentially on continental plates
    let mut rng_state = ov.simple_seed;
    
    for continent_id in 0..continents_n {
        // Pick a continental plate for this continent
        if !continental_plates.is_empty() {
            rng_state = rng_state.wrapping_mul(1664525).wrapping_add(1013904223);
            let plate_idx = (rng_state % (continental_plates.len() as u64)) as usize;
            let target_plate_id = continental_plates[plate_idx];
            
            // Find a good center on this continental plate
            let center_idx = find_continental_plate_center(world, target_plate_id, rng_state);
            
            // Generate continent size based on preset
            let continent_area_fraction = match ov.simple_preset {
                1 => 0.15, // Archipelago - smaller continents
                2 => 0.4,  // Continental - larger continents  
                3 => 0.8,  // Supercontinent - one massive continent
                _ => 0.25, // Balanced
            };
            
            let continent_area = target_land_area_m2 * continent_area_fraction / (continents_n as f64);
            
            // Generate realistic continent with multifractal noise
            generate_continent_with_noise(world, center_idx, continent_area, continent_id, rng_state, target_plate_id);
            
            println!("[world_gen] Continent {} placed on plate {} at cell {}", continent_id, target_plate_id, center_idx);
        }
    }
}

/// Classify plates as continental or oceanic based on size and position
fn classify_continental_plates(world: &engine::world::World, _ov: &overlay::OverlayState) -> Vec<u16> {
    let mut plate_areas = std::collections::HashMap::new();
    let mut plate_centers = std::collections::HashMap::new();
    let mut plate_counts = std::collections::HashMap::new();
    
    // Calculate area and center for each plate
    for i in 0..world.grid.cells {
        let plate_id = world.plates.plate_id[i];
        let area = world.area_m2[i] as f64;
        let pos = world.grid.pos_xyz[i];
        
        *plate_areas.entry(plate_id).or_insert(0.0) += area;
        *plate_counts.entry(plate_id).or_insert(0) += 1;
        
        let center = plate_centers.entry(plate_id).or_insert([0.0f64; 3]);
        center[0] += pos[0] as f64;
        center[1] += pos[1] as f64;
        center[2] += pos[2] as f64;
    }
    
    // Normalize centers and classify plates
    let mut continental_plates = Vec::new();
    let total_area: f64 = world.area_m2.iter().map(|&a| a as f64).sum();
    let mean_plate_area = total_area / (plate_areas.len() as f64);
    
    for (&plate_id, &area) in &plate_areas {
        // Normalize center
        if let Some(center) = plate_centers.get_mut(&plate_id) {
            let count = plate_counts[&plate_id] as f64;
            center[0] /= count;
            center[1] /= count; 
            center[2] /= count;
        }
        
            // Classify as continental if larger than average (continental plates are typically larger)
        if area > mean_plate_area * 0.8 {
            continental_plates.push(plate_id);
            println!("[world_gen] Plate {} classified as CONTINENTAL (area={:.1e} m², mean={:.1e})", 
                plate_id, area, mean_plate_area);
        } else {
            println!("[world_gen] Plate {} classified as oceanic (area={:.1e} m², mean={:.1e})", 
                plate_id, area, mean_plate_area);
        }
    }
    
    // Ensure we have at least some continental plates
    if continental_plates.is_empty() {
        // If no large plates, pick the largest few plates
        let mut plate_list: Vec<_> = plate_areas.iter().collect();
        plate_list.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        
        let num_continental = (plate_areas.len() / 2).max(1).min(4); // 50% of plates, 1-4 range
        for i in 0..num_continental {
            if i < plate_list.len() {
                continental_plates.push(*plate_list[i].0);
            }
        }
    }
    
    continental_plates
}

/// Find a good center point on a specific continental plate
fn find_continental_plate_center(world: &engine::world::World, target_plate_id: u16, mut rng_state: u64) -> usize {
    // Collect all cells belonging to this plate
    let mut plate_cells = Vec::new();
    for i in 0..world.grid.cells {
        if world.plates.plate_id[i] == target_plate_id {
            plate_cells.push(i);
        }
    }
    
    if plate_cells.is_empty() {
        return 0; // Fallback
    }
    
    // Pick a random cell from this plate
    rng_state = rng_state.wrapping_mul(1664525).wrapping_add(1013904223);
    let idx = (rng_state % (plate_cells.len() as u64)) as usize;
    plate_cells[idx]
}

/// Generate a single continent with realistic multifractal noise
fn generate_continent_with_noise(
    world: &mut engine::world::World,
    center_idx: usize,
    target_area: f64,
    continent_id: u32,
    rng_state: u64,
    target_plate_id: u16,
) {
    if center_idx >= world.grid.cells {
        return;
    }
    
    let center_pos = world.grid.pos_xyz[center_idx];
    
    // Estimate continent radius from target area (spherical approximation)
    let mean_cell_area = world.area_m2.iter().map(|&a| a as f64).sum::<f64>() / (world.grid.cells as f64);
    let target_cells = target_area / mean_cell_area;
    let continent_radius = (target_cells / std::f64::consts::PI).sqrt() * 0.3; // Increased scale for visibility
    
    println!("[world_gen] Continent {}: target_area={:.1e} m², target_cells={:.0}, radius={:.3}", 
        continent_id, target_area, target_cells, continent_radius);
    
    // Multifractal noise parameters
    let noise_scales = [0.5, 1.0, 2.0, 4.0, 8.0]; // Multiple octaves
    let noise_weights = [1.0, 0.5, 0.25, 0.125, 0.0625]; // Decreasing weights
    
    for i in 0..world.grid.cells {
        // Only place continents on the target continental plate
        if world.plates.plate_id[i] != target_plate_id {
            continue;
        }
        
        let pos = world.grid.pos_xyz[i];
        let distance = ((pos[0] - center_pos[0]).powi(2) + 
                       (pos[1] - center_pos[1]).powi(2) + 
                       (pos[2] - center_pos[2]).powi(2)).sqrt();
        
        if distance < (continent_radius * 2.0) as f32 {
            // Base falloff from center
            let base_falloff = (1.0 - (distance as f64 / (continent_radius * 1.5))).max(0.0);
            
            if base_falloff > 0.0 {
                // Generate multifractal noise for realistic continental shape
                let mut noise_value = 0.0;
                let mut total_weight = 0.0;
                
                for (scale, weight) in noise_scales.iter().zip(noise_weights.iter()) {
                    let noise_freq = scale / continent_radius;
                    let noise_sample = simplex_noise_3d(
                        pos[0] as f64 * noise_freq,
                        pos[1] as f64 * noise_freq,
                        pos[2] as f64 * noise_freq,
                        rng_state.wrapping_add(continent_id as u64),
                    );
                    noise_value += noise_sample * weight;
                    total_weight += weight;
                }
                
                noise_value /= total_weight;
                
                // Combine base falloff with noise for realistic continental margins
                let continental_strength = (base_falloff * (0.7 + 0.3 * noise_value)).max(0.0);
                
                if continental_strength > 0.05 { // Lower threshold to ensure continents are created
                    // Continental fraction (0.6-0.9 based on noise)
                    let c_value = (0.6 + 0.3 * continental_strength) as f32;
                    world.c[i] = c_value.max(world.c[i]);
                    
                    // Continental thickness variation (20-45 km based on noise and position)
                    let thickness_base = 30_000.0; // 30 km base (more realistic)
                    let thickness_variation = 20_000.0 * continental_strength; // Up to 20 km variation
                    let thickness = (thickness_base + thickness_variation) as f32;
                    world.th_c_m[i] = thickness.max(world.th_c_m[i]);
                    
                    // Debug: Log first few continental cells
                    if continent_id == 0 && world.c[i] > 0.5 {
                        println!("[world_gen] Continental cell {}: C={:.2}, th_c={:.0}m, strength={:.2}", 
                                i, world.c[i], world.th_c_m[i], continental_strength);
                    }
                }
            }
        }
    }
}

/// Add realistic topographic variation using multifractal noise
fn add_realistic_topography(world: &mut engine::world::World, ov: &overlay::OverlayState) {
    let rng_state = ov.simple_seed.wrapping_add(12345); // Different seed for topography
    
    for i in 0..world.grid.cells {
        let pos = world.grid.pos_xyz[i];
        
        // Only add topography to continental areas
        if world.c[i] > 0.1 {
            // Generate multifractal elevation noise
            let mut elevation_noise = 0.0;
            let mut total_weight = 0.0;
            
            // Multiple octaves for realistic mountain ranges
            let scales = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
            let weights = [1.0, 0.6, 0.4, 0.3, 0.2, 0.1];
            
            for (scale, weight) in scales.iter().zip(weights.iter()) {
                let noise_freq = scale * 0.1; // Scale for continental features
                let noise_sample = simplex_noise_3d(
                    pos[0] as f64 * noise_freq,
                    pos[1] as f64 * noise_freq, 
                    pos[2] as f64 * noise_freq,
                    rng_state,
                );
                elevation_noise += noise_sample * weight;
                total_weight += weight;
            }
            
            elevation_noise /= total_weight;
            
            // Scale elevation based on continental strength and thickness
            let continental_strength = world.c[i];
            let max_elevation = 3000.0; // 3 km maximum elevation
            let elevation_variation = (elevation_noise * max_elevation * continental_strength as f64) as f32;
            
            // Apply topographic variation to depth (negative for elevation)
            world.depth_m[i] -= elevation_variation.abs(); // Make sure land is above sea level
        }
    }
}

/// Simplified 3D simplex noise implementation
fn simplex_noise_3d(x: f64, y: f64, z: f64, seed: u64) -> f64 {
    // Simple deterministic pseudo-noise using trigonometric functions
    // This is a simplified version - in production you'd use a proper noise library
    
    let mut hash = seed;
    hash = hash.wrapping_mul(1664525).wrapping_add(1013904223);
    
    let a = (x * 12.9898 + y * 78.233 + z * 37.719 + hash as f64 * 0.001).sin() * 43758.5453;
    let b = (x * 93.9898 + y * 67.345 + z * 83.217 + hash as f64 * 0.002).sin() * 47896.6789;
    let c = (x * 23.1234 + y * 56.789 + z * 91.234 + hash as f64 * 0.003).sin() * 39217.9876;
    
    // Combine and normalize to [-1, 1]
    let noise = (a.fract() + b.fract() + c.fract()) / 3.0;
    (noise - 0.5) * 2.0
}

/// Simple fallback continent generation that definitely works
fn create_simple_fallback_continents(
    world: &mut engine::world::World,
    ov: &overlay::OverlayState,
    continents_n: u32,
) {
    println!("[world_gen] Creating {} fallback continents", continents_n);
    
    let n_cells = world.grid.cells;
    let mut rng_state = ov.simple_seed;
    
    // Create simple, guaranteed-to-work continents
    for continent_id in 0..continents_n {
        rng_state = rng_state.wrapping_mul(1664525).wrapping_add(1013904223);
        let center_idx = (rng_state % (n_cells as u64)) as usize;
        
        // Simple circular continent - guaranteed to create material
        let center_pos = world.grid.pos_xyz[center_idx];
        let radius = 0.2f32; // Large enough to be visible
        
        let mut cells_modified = 0;
        for i in 0..n_cells {
            let pos = world.grid.pos_xyz[i];
            let distance = ((pos[0] - center_pos[0]).powi(2) + 
                           (pos[1] - center_pos[1]).powi(2) + 
                           (pos[2] - center_pos[2]).powi(2)).sqrt();
            
            if distance < radius {
                let falloff = (1.0 - (distance / radius)).max(0.0);
                if falloff > 0.1 {
                    world.c[i] = (0.8 * falloff).max(world.c[i]);
                    world.th_c_m[i] = (40_000.0 * falloff).max(world.th_c_m[i]); // 40 km thick
                    cells_modified += 1;
                }
            }
        }
        
        println!("[world_gen] Fallback continent {} at cell {}: radius={:.2}, cells_modified={}", 
            continent_id, center_idx, radius, cells_modified);
    }
    
    let c_total: f32 = world.c.iter().sum();
    let th_total: f32 = world.th_c_m.iter().sum();
    println!("[world_gen] Fallback result: C_sum={:.1}, th_c_sum={:.1} m", c_total, th_total);
}

/// 🎨 Visual Overlays - Organized by complexity with keyboard shortcuts
fn render_visual_overlays_unified(
    ui: &mut egui::Ui,
    _ctx: &egui::Context,
    _world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
) {
    // ============= FLAT UI STRUCTURE - NO NESTED GROUPS =============
    
    // Visual Overlays Section
    egui::CollapsingHeader::new("🎨 Visual Overlays")
        .default_open(false)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.show_plates, "Plate Boundaries [1]");
                ui.checkbox(&mut ov.show_vel, "Plate Velocity [2]");
                ui.checkbox(&mut ov.show_continents, "Continental Crust [C]");
                ui.checkbox(&mut ov.show_bounds, "Boundary Types [3]");
            });
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.show_subduction, "Subduction Zones [7]");
                ui.checkbox(&mut ov.show_transforms, "Transform Faults [0]");
                ui.checkbox(&mut ov.show_age_depth, "Age-Depth [6]");
                ui.checkbox(&mut ov.show_plate_id, "Plate ID Colors");
            });
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.show_triple_junctions, "Triple Junctions [5]");
                ui.checkbox(&mut ov.show_plate_adjacency, "Plate Adjacency [4]");
                ui.checkbox(&mut ov.show_flexure, "Flexure Response");
                ui.checkbox(&mut ov.debug_wireframes, "Debug Wireframes");
            });
            
            ui.horizontal(|ui| {
                ui.label("Color mode:");
                ui.selectable_value(&mut ov.color_mode, 0, "Elevation");
                ui.selectable_value(&mut ov.color_mode, 1, "Biomes");
                ui.checkbox(&mut ov.shade_on, "Hillshading");
                if ov.shade_on {
                    ui.add(egui::Slider::new(&mut ov.shade_strength, 0.0..=1.0).text("Strength"));
                }
                ui.checkbox(&mut ov.legend_on, "Show Legend");
            });
            
            // Invalidate caches when overlays change
            if ui.ui_contains_pointer() {
                ov.plates_cache = None;
                ov.vel_cache = None;
                ov.bounds_cache = None;
                ov.color_dirty = true;
            }
        });

    // Force Balance Parameters Section
    egui::CollapsingHeader::new("⚖️ Force Balance Parameters")
        .default_open(false)
        .show(ui, |ui| {
            ui.checkbox(&mut ov.enable_force_balance, "Enable Force Balance")
                .on_hover_text("Controls automatic adjustment of plate rotation rates based on boundary forces");
            
            ui.label("Gain (dimensionless):");
            ui.add(egui::Slider::new(&mut ov.fb_gain, 1e-12..=1e-5).logarithmic(true))
                .on_hover_text("Primary control for plate motion responsiveness. Higher = faster rotation changes. Typical: 1e-10 to 1e-8");
            
            ui.label("Damping (per Myr):");
            ui.add(egui::Slider::new(&mut ov.fb_damp_per_myr, 0.001..=1.0).logarithmic(true))
                .on_hover_text("Stabilizes plate motion. Higher = more damping, slower changes. Typical: 0.01 to 0.1");
            
            ui.label("Convergent Coefficient:");
            ui.add(egui::Slider::new(&mut ov.fb_k_conv, 0.1..=5.0))
                .on_hover_text("Force multiplier for convergent boundaries (subduction/collision). Typical: 0.5 to 2.0");
            
            ui.label("Divergent Coefficient:");
            ui.add(egui::Slider::new(&mut ov.fb_k_div, 0.1..=5.0))
                .on_hover_text("Force multiplier for divergent boundaries (spreading ridges). Typical: 0.2 to 1.0");
            
            ui.label("Transform Coefficient:");
            ui.add(egui::Slider::new(&mut ov.fb_k_trans, 0.01..=1.0))
                .on_hover_text("Force multiplier for transform boundaries (lateral motion). Typical: 0.05 to 0.2");
            
            ui.label("Max Ω Change (rad/yr per step):");
            ui.add(egui::Slider::new(&mut ov.fb_max_domega, 1e-12..=1e-6).logarithmic(true))
                .on_hover_text("Maximum rotation rate change per time step. Prevents instability. Typical: 1e-9 to 1e-7");
            
            ui.label("Max Ω Total (rad/yr):");
            ui.add(egui::Slider::new(&mut ov.fb_max_omega, 1e-9..=1e-4).logarithmic(true))
                .on_hover_text("Maximum absolute plate rotation rate. Real plates: ~1e-7 to 1e-5 rad/yr");
            
            ui.label("Force Balance Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_force_balance_every, 1..=10))
                .on_hover_text("How often to apply force balance updates. 1 = every step (dynamic), higher = less frequent. CRITICAL: Must be >0 for plates to move!");
        });

    // Plate Motion & Kinematics
    egui::CollapsingHeader::new("🔄 Plate Motion & Kinematics")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Rigid Motion Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_rigid_motion, 1..=10))
                .on_hover_text("How often to update plate velocities. 1=every step, higher=less frequent. Affects motion smoothness.");
            
            ui.label("Force Balance Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_force_balance, 1..=20))
                .on_hover_text("How often to adjust plate rotation rates. Higher = more stable but less responsive motion.");
        });

    // Plate Boundaries & Tectonics  
    egui::CollapsingHeader::new("🌋 Plate Boundaries & Tectonics")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Transform Faults Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_transforms, 1..=20))
                .on_hover_text("Frequency of transform fault processing. Controls lateral plate motion and strike-slip zones.");
            
            ui.label("Subduction Zones Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_subduction, 1..=20))
                .on_hover_text("Frequency of subduction processing. Controls convergent margin dynamics and volcanic arcs.");
            
            ui.label("Oceanic Rifting Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_rifting, 0..=50))
                .on_hover_text("Continental breakup and oceanic rifting. 0=disabled. Controls new ocean formation.");
            
            ui.label("Ridge Birth Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_ridge_birth, 0..=50))
                .on_hover_text("New spreading ridge formation. 0=disabled. Resets seafloor age at divergent boundaries.");
        });

    // Mountain Building & Collisions
    egui::CollapsingHeader::new("🏔️ Mountain Building & Collisions")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Orogeny Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_orogeny, 1..=50))
                .on_hover_text("Continental collision and mountain building. Creates orogenic belts when continents collide.");
            
            ui.label("Oceanic Accretion Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_accretion, 0..=50))
                .on_hover_text("Growth of volcanic arcs and accretionary wedges at subduction zones. 0=disabled.");
        });

    // Isostasy & Buoyancy
    egui::CollapsingHeader::new("⚖️ Isostasy & Buoyancy")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Flexural Response Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_flexure, 1..=10))
                .on_hover_text("Lithospheric bending under loads. Controls crustal deflection from mountains and sediments.");
            
            ui.label("Isostatic Adjustment Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_isostasy, 1..=10))
                .on_hover_text("Global sea level regulation and isostatic rebound. Maintains target land fraction.");
            
            ui.label("Continental Buoyancy Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_continental_buoyancy, 1..=10))
                .on_hover_text("Buoyant response of continental crust. Controls continental elevation relative to oceans.");
        });

    // Surface Processes
    egui::CollapsingHeader::new("🌊 Surface Processes")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Surface Processes Cadence (steps):");
            ui.add(egui::Slider::new(&mut ov.cadence_surface_processes, 1..=50))
                .on_hover_text("Erosion, sedimentation, and landscape evolution. Controls weathering and transport of material.");
        });
}

/// ⚙️ Physics Processes - Progressive complexity with smart defaults
fn render_physics_processes_unified(
    ui: &mut egui::Ui,
    _ctx: &egui::Context,
    _world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
) {
    egui::CollapsingHeader::new("⚙️ Physics Processes")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Core Processes (Essential):");
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.enable_rigid_motion, "Rigid Motion")
                    .on_hover_text("✅ Essential: plate movement and kinematics");
                ui.checkbox(&mut ov.enable_flexure, "Flexure")
                    .on_hover_text("✅ Essential: crustal bending under loads");
            });
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.enable_isostasy, "Isostasy")
                    .on_hover_text("✅ Essential: sea level regulation and crustal equilibrium");
                ui.checkbox(&mut ov.enable_continental_buoyancy, "Continental Buoyancy")
                    .on_hover_text("✅ Essential: continent elevation above sea level");
            });
            
            ui.separator();
            
            ui.label("Geological Processes:");
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.enable_subduction, "Subduction")
                    .on_hover_text("Convergent plate boundaries - where oceanic plates dive under others");
                ui.checkbox(&mut ov.enable_transforms, "Transforms")
                    .on_hover_text("Lateral plate motion - strike-slip faults and shear zones");
            });
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.enable_surface_processes, "Surface Processes")
                    .on_hover_text("Erosion and sediment transport - shapes landscapes over time");
                ui.checkbox(&mut ov.enable_orogeny, "Orogeny")
                    .on_hover_text("Mountain building from continent-continent collision");
            });
            
            // Advanced Processes (collapsible)
            ui.collapsing("🔬 Advanced Processes", |ui| {
                ui.horizontal(|ui| {
                    ui.checkbox(&mut ov.enable_accretion, "Accretion")
                        .on_hover_text("Expert: terrane attachment at convergent margins");
                    ui.checkbox(&mut ov.enable_rifting, "Rifting")
                        .on_hover_text("Expert: continental breakup and passive margin formation");
                });
                
                ui.horizontal(|ui| {
                    ui.checkbox(&mut ov.enable_ridge_birth, "Ridge Birth")
                        .on_hover_text("Expert: formation of new spreading centers");
                    ui.checkbox(&mut ov.enable_force_balance, "Force Balance")
                        .on_hover_text("Expert: plate motion dynamics and driving forces");
                });
            });
            
            ui.separator();
            
            // Backend Configuration (collapsible)
            ui.collapsing("🔬 Backend Configuration", |ui| {
                ui.horizontal(|ui| {
                    ui.label("Flexure:")
                        .on_hover_text("Choose computation method for crustal flexure");
                    if ui.selectable_label(ov.flexure_backend_cpu, "CPU (Winkler)")
                        .on_hover_text("✅ Reliable - Always works, good for learning")
                        .clicked() {
                        ov.flexure_backend_cpu = true;
                    }
                    if ui.selectable_label(!ov.flexure_backend_cpu, "GPU (Multigrid)")
                        .on_hover_text("🔬 Advanced - Faster but requires compatible hardware")
                        .clicked() {
                        ov.flexure_backend_cpu = false;
                    }
                });
                
                if !ov.flexure_backend_cpu {
                    ui.horizontal(|ui| {
                        ui.label("GPU Levels:");
                        ui.add(egui::Slider::new(&mut ov.flexure_gpu_levels, 1..=5)
                            .show_value(true))
                            .on_hover_text("Multigrid levels - more levels = more accurate but slower");
                        ui.label("V-cycles:");
                        ui.add(egui::Slider::new(&mut ov.flexure_gpu_cycles, 1..=8)
                            .show_value(true))
                            .on_hover_text("Solver iterations - more cycles = more accurate");
                    });
                }
            });
            
            // Performance Tuning (collapsible)
            ui.collapsing("🔬 Performance Tuning", |ui| {
                ui.horizontal(|ui| {
                    ui.label("Cadence preset:")
                        .on_hover_text("Balance between simulation quality and performance");
                    ui.selectable_value(&mut ov.cadence_preset, "Balanced".to_string(), "Balanced")
                        .on_hover_text("✅ Recommended - good quality and performance");
                    ui.selectable_value(&mut ov.cadence_preset, "Performance".to_string(), "Performance")
                        .on_hover_text("⚡ Faster - reduced quality for speed");
                    ui.selectable_value(&mut ov.cadence_preset, "Quality".to_string(), "Quality")
                        .on_hover_text("🔬 Slower - maximum quality and accuracy");
                });
                
                if ov.cadence_preset == "Custom" {
                    ui.label("Custom cadences (steps between executions):");
                    // Custom cadence controls would go here
                    ui.label("(Custom cadence controls - see Process Cadences panel)");
                }
            });
        });
}

/// 🐛 Debug & Diagnostics - Expert features and developer tools
fn render_debug_diagnostics_unified(
    ui: &mut egui::Ui,
    _ctx: &egui::Context,
    _world: &mut engine::world::World,
    ov: &mut overlay::OverlayState,
) {
    egui::CollapsingHeader::new("🐛 Debug & Diagnostics")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("🔬 Expert Features:");
            
            ui.horizontal(|ui| {
                ui.checkbox(&mut ov.debug_enable_all, "Run All Passes")
                    .on_hover_text("Force all physics processes to run every step (ignores cadences)");
                ui.checkbox(&mut ov.disable_subduction, "Disable Subduction")
                    .on_hover_text("Temporarily disable subduction for debugging");
            });
            
            ui.separator();
            
            ui.label("Performance & Logging:");
            ui.horizontal(|ui| {
                ui.label("Profiling burst:");
                ui.add(egui::DragValue::new(&mut ov.debug_burst_steps).speed(1))
                    .on_hover_text("Number of steps to profile for performance analysis");
            });
            
            ui.separator();
            
            ui.label("Export Tools:");
            if ui.button("📄 Export Debug CSV")
                .on_hover_text("Export current simulation state for analysis")
                .clicked() {
                // Export functionality would be triggered here
            }
            
            if ui.button("📊 Export Performance Log")
                .on_hover_text("Export timing and performance metrics")
                .clicked() {
                // Performance export would be triggered here  
            }
        });
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
        // Prefer DX12 on Windows to avoid Vulkan present-mode spam; allow env override
        let mut backends = wgpu::Backends::DX12;
        if std::env::var_os("WGPU_BACKEND").is_some() {
            backends = wgpu::util::backend_bits_from_env().unwrap_or(wgpu::Backends::DX12);
        }
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends,
            dx12_shader_compiler: wgpu::Dx12Compiler::Fxc,
            flags: wgpu::InstanceFlags::DEBUG,
            gles_minor_version: wgpu::Gles3MinorVersion::Automatic,
        });
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
            .find(|f| {
                matches!(f, wgpu::TextureFormat::Bgra8Unorm | wgpu::TextureFormat::Rgba8Unorm)
            })
            .unwrap_or_else(|| {
                surface_caps
                    .formats
                    .iter()
                    .copied()
                    .find(|f| !f.is_srgb())
                    .unwrap_or(surface_caps.formats[0])
            });

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
    // Initialize tracing subscriber with env filter; default to info if not set
    {
        use tracing_subscriber::{fmt, EnvFilter};
        let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
        let _ = fmt::Subscriber::builder().with_env_filter(filter).try_init();
    }
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
    // GPU buffer manager for consistent elevation data
    let mut gpu_buf_mgr = GpuBufferManager::new();
    // Edge-triggered snapshots state
    let _next_snapshot_t: f64 = f64::INFINITY;
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

    // Background simulation thread with unified world state sync
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{mpsc, Arc, Mutex};
    
    // Single world update channel - no separate elevation stream
    let (tx_world, rx_world) = mpsc::channel::<WorldSnapshot>();
    let (tx_cmd, rx_cmd) = mpsc::channel::<SimCommand>();
    let sim_busy = Arc::new(AtomicBool::new(false));
    let sim_handle: Arc<Mutex<Option<std::thread::JoinHandle<()>>>> = Arc::new(Mutex::new(None));
    
    // Initialize background simulation thread
    {
        let tx_world_bg = tx_world.clone();
        let sim_busy_bg = sim_busy.clone();
        let rx_cmd_bg = rx_cmd;
        let init_world = WorldSnapshot {
            depth_m: world.depth_m.clone(),
            c: world.c.clone(),
            th_c_m: world.th_c_m.clone(),
            sea_eta_m: world.sea.eta_m,
            plate_id: world.plates.plate_id.clone(),
            pole_axis: world.plates.pole_axis.clone(),
            omega_rad_yr: world.plates.omega_rad_yr.clone(),
            v_en: world.v_en.clone(),
            age_myr: world.age_myr.clone(),
        };
        
        let sim_thread = std::thread::spawn(move || {
            let mut sim_world = engine::world::World::new(64, 8, 12345);
            
            // Initialize from snapshot
            if sim_world.depth_m.len() == init_world.depth_m.len() {
                sim_world.depth_m = init_world.depth_m;
                sim_world.depth_stage_m.clone_from(&sim_world.depth_m);
            }
            if sim_world.c.len() == init_world.c.len() { sim_world.c = init_world.c; }
            if sim_world.th_c_m.len() == init_world.th_c_m.len() { sim_world.th_c_m = init_world.th_c_m; }
            sim_world.sea.eta_m = init_world.sea_eta_m;
            if sim_world.plates.plate_id.len() == init_world.plate_id.len() { sim_world.plates.plate_id = init_world.plate_id; }
            if sim_world.plates.pole_axis.len() == init_world.pole_axis.len() { sim_world.plates.pole_axis = init_world.pole_axis; }
            if sim_world.plates.omega_rad_yr.len() == init_world.omega_rad_yr.len() { sim_world.plates.omega_rad_yr = init_world.omega_rad_yr; }
            if sim_world.v_en.len() == init_world.v_en.len() { sim_world.v_en = init_world.v_en; }
            if sim_world.age_myr.len() == init_world.age_myr.len() { sim_world.age_myr = init_world.age_myr; }
            sim_world.boundaries = engine::boundaries::Boundaries::classify(&sim_world.grid, &sim_world.plates.plate_id, &sim_world.v_en, 0.005);
            
            while let Ok(cmd) = rx_cmd_bg.recv() {
                match cmd {
                    SimCommand::Step(cfg, process_flags) => {
                        sim_busy_bg.store(true, Ordering::SeqCst);
                        
                        // Convert UI flags to unified PhysicsConfig
                        let mut config = engine::config::PhysicsConfig::simple_mode();
                        config.dt_myr = cfg.dt_myr;
                        config.target_land_frac = cfg.target_land_frac;
                        config.freeze_eta = cfg.freeze_eta;
                        
                        // Apply process enables
                        config.enable_rigid_motion = process_flags.enable_rigid_motion;
                        config.enable_subduction = process_flags.enable_subduction;
                        config.enable_transforms = process_flags.enable_transforms;
                        config.enable_flexure = process_flags.enable_flexure;
                        config.enable_surface_processes = process_flags.enable_surface_processes;
                        config.enable_isostasy = process_flags.enable_isostasy;
                        config.enable_continental_buoyancy = process_flags.enable_continental_buoyancy;
                        config.enable_orogeny = process_flags.enable_orogeny;
                        config.enable_accretion = process_flags.enable_accretion;
                        config.enable_rifting = process_flags.enable_rifting;
                        config.enable_ridge_birth = process_flags.enable_ridge_birth;
                        
                        // Apply flexure backend
                        config.flexure_backend = if process_flags.flexure_backend_cpu {
                            engine::flexure_manager::FlexureBackend::CpuWinkler
                        } else {
                            engine::flexure_manager::FlexureBackend::GpuMultigrid {
                                levels: process_flags.flexure_gpu_levels,
                                cycles: process_flags.flexure_gpu_cycles,
                            }
                        };
                        
                        // Apply cadence config
                        config.cadence_config = process_flags.cadence_config;
                        
                        // Run unified pipeline step
                        let mut pipeline = engine::unified_pipeline::UnifiedPipeline::new(config);
                        let mode = engine::config::PipelineMode::Realtime { preserve_depth: true };
                        let _result = pipeline.step(&mut sim_world, mode);
                        
                        // Send complete world state (unified approach)
                        let ws = WorldSnapshot {
                            depth_m: sim_world.depth_m.clone(),
                            c: sim_world.c.clone(),
                            th_c_m: sim_world.th_c_m.clone(),
                            sea_eta_m: sim_world.sea.eta_m,
                            plate_id: sim_world.plates.plate_id.clone(),
                            pole_axis: sim_world.plates.pole_axis.clone(),
                            omega_rad_yr: sim_world.plates.omega_rad_yr.clone(),
                            v_en: sim_world.v_en.clone(),
                            age_myr: sim_world.age_myr.clone(),
                        };
                        let _ = tx_world_bg.send(ws);
                        sim_busy_bg.store(false, Ordering::SeqCst);
                    }
                    SimCommand::SyncWorld(ws) => {
                        // Update simulation world from UI (e.g., after Generate World)
                        if sim_world.depth_m.len() == ws.depth_m.len() {
                            sim_world.depth_m = ws.depth_m;
                            sim_world.depth_stage_m.clone_from(&sim_world.depth_m);
                        }
                        if sim_world.c.len() == ws.c.len() { sim_world.c = ws.c; }
                        if sim_world.th_c_m.len() == ws.th_c_m.len() { sim_world.th_c_m = ws.th_c_m; }
                        sim_world.sea.eta_m = ws.sea_eta_m;
                        if sim_world.plates.plate_id.len() == ws.plate_id.len() { sim_world.plates.plate_id = ws.plate_id; }
                        if sim_world.plates.pole_axis.len() == ws.pole_axis.len() { sim_world.plates.pole_axis = ws.pole_axis; }
                        if sim_world.plates.omega_rad_yr.len() == ws.omega_rad_yr.len() { sim_world.plates.omega_rad_yr = ws.omega_rad_yr; }
                        if sim_world.v_en.len() == ws.v_en.len() { sim_world.v_en = ws.v_en; } else {
                            sim_world.v_en = engine::plates::velocity_en_m_per_yr(&sim_world.grid, &sim_world.plates, &sim_world.plates.plate_id);
                        }
                        if sim_world.age_myr.len() == ws.age_myr.len() { sim_world.age_myr = ws.age_myr; }
                        sim_world.boundaries = engine::boundaries::Boundaries::classify(&sim_world.grid, &sim_world.plates.plate_id, &sim_world.v_en, 0.005);
                    }
                    SimCommand::Stop => break,
                }
            }
        });
        
        if let Ok(mut h) = sim_handle.lock() {
            *h = Some(sim_thread);
        }
    }

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
                    WindowEvent::CloseRequested => {
                        // Graceful shutdown of simulation thread
                        let _ = tx_cmd.send(SimCommand::Stop);
                        if let Ok(mut h) = sim_handle.lock() {
                            if let Some(jh) = h.take() {
                                let _ = jh.join();
                            }
                        }
                        elwt.exit()
                    },
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
                            // Unified overlay shortcuts available in both Simple and Advanced modes
                                if ctx.input(|i| i.key_pressed(egui::Key::Num1)) { ov.show_plates = !ov.show_plates; ov.plates_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num2)) { ov.show_vel = !ov.show_vel; ov.vel_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num3)) { ov.show_bounds = !ov.show_bounds; ov.bounds_cache = None; }
                            if ctx.input(|i| i.key_pressed(egui::Key::C)) { ov.show_continents = !ov.show_continents; if ov.show_continents && (ov.mesh_continents.is_none() || ov.mesh_coastline.is_none()) { _continents_dirty = true; } }
                            if ctx.input(|i| i.key_pressed(egui::Key::H)) { ov.show_hud = !ov.show_hud; }
                            
                            // Advanced-only shortcuts (complex overlays and debug features)
                            if ov.ui_mode.show_advanced_overlays() {
                                if ctx.input(|i| i.key_pressed(egui::Key::Num4)) { ov.show_plate_adjacency = !ov.show_plate_adjacency; ov.net_adj_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num5)) { ov.show_triple_junctions = !ov.show_triple_junctions; ov.net_tj_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num6)) { ov.show_age_depth = !ov.show_age_depth; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num7)) { ov.show_subduction = !ov.show_subduction; ov.subd_trench=None; ov.subd_arc=None; ov.subd_backarc=None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Num0)) { ov.show_transforms = !ov.show_transforms; ov.trans_pull=None; ov.trans_rest=None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::A)) { age_plot.show = !age_plot.show; }
                                if ctx.input(|i| i.key_pressed(egui::Key::M)) { ov.show_map_color_panel = !ov.show_map_color_panel; }
                                if ctx.input(|i| i.key_pressed(egui::Key::S)) { ov.surface_enable = !ov.surface_enable; }
                                if ctx.input(|i| i.key_pressed(egui::Key::L)) { ov.apply_sea_level = !ov.apply_sea_level; ov.bathy_cache = None; }
                                if ctx.input(|i| i.key_pressed(egui::Key::F)) { flex.show = !flex.show; if flex.show { flex.recompute(); } }
                                if ctx.input(|i| i.key_pressed(egui::Key::G)) { ov.show_flexure = !ov.show_flexure; _flex_dirty = true; }
                                if ctx.input(|i| i.key_pressed(egui::Key::Y)) { ov.show_hypsometry = !ov.show_hypsometry; }
                            }
                            
                            // GPU raster provides base elevation colors; disable redundant CPU bathy overlay
                            ov.show_bathy = !ov.use_gpu_raster;

                            egui::TopBottomPanel::top("hud").show(ctx, |ui| {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label(format!("Aulé Viewer v{}", env!("CARGO_PKG_VERSION")));
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
                                    // Live land fraction (area-weighted) directly from world data
                                    let land_pct_now: f64 = {
                                        let mut area_sum = 0.0f64; 
                                        let mut land_area = 0.0f64;
                                        for (&d, &a) in world.depth_m.iter().zip(world.area_m2.iter()) { 
                                            area_sum += a as f64; 
                                            if (world.sea.eta_m as f64 - d as f64) > 0.0 { land_area += a as f64; } 
                                        }
                                        if area_sum > 0.0 { 100.0 * (land_area / area_sum) } else { 0.0 }
                                    };
                                    ui.label(format!("Land: {:.1}%", land_pct_now));
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
                                                // Render unified progressive UI (replaces Simple/Advanced modes)
                                                render_unified_progressive_panels(ui, ctx, &mut world, &mut ov, &tx_cmd);
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
                            // Atomic world state updates - drain ALL pending updates before rendering
                            let mut latest_world_update: Option<WorldSnapshot> = None;
                            while let Ok(ws) = rx_world.try_recv() {
                                latest_world_update = Some(ws); // Keep only the most recent update
                            }
                            
                            // Apply the most recent world state atomically (if any)
                            if let Some(ws) = latest_world_update {
                                // Update complete world state from simulation thread
                                if world.depth_m.len() == ws.depth_m.len() { world.depth_m = ws.depth_m; }
                                if world.c.len() == ws.c.len() { world.c = ws.c; }
                                if world.th_c_m.len() == ws.th_c_m.len() { world.th_c_m = ws.th_c_m; }
                                world.sea.eta_m = ws.sea_eta_m;
                                world.clock.t_myr = ws.age_myr.iter().fold(0.0, |acc, &age| acc.max(age as f64)); // Approximate time from max age
                                if world.plates.plate_id.len() == ws.plate_id.len() { world.plates.plate_id = ws.plate_id; }
                                if world.plates.pole_axis.len() == ws.pole_axis.len() { world.plates.pole_axis = ws.pole_axis; }
                                if world.plates.omega_rad_yr.len() == ws.omega_rad_yr.len() { world.plates.omega_rad_yr = ws.omega_rad_yr; }
                                if world.v_en.len() == ws.v_en.len() { world.v_en = ws.v_en; }
                                if world.age_myr.len() == ws.age_myr.len() { world.age_myr = ws.age_myr; }
                                
                                // Recompute boundaries with updated state
                                world.boundaries = engine::boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
                                
                                // Mark everything dirty for immediate update
                                ov.world_dirty = true; ov.color_dirty = true;
                                ov.plates_cache = None;
                                ov.vel_cache = None;
                                ov.bounds_cache = None;
                                ov.age_cache = None;
                                ov.bathy_cache = None;
                                ov.net_adj_cache = None;
                                ov.net_tj_cache = None;
                                // Invalidate GPU buffer cache when world changes
                                gpu_buf_mgr.invalidate_cache();
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
                                        // Upload initial heights and LUT (render-only clamped copy)
                                        let heights_init: Vec<f32> = elevation_curr_clone()
                                            .unwrap_or_else(|| world.depth_m.iter().map(|&d| world.sea.eta_m - d).collect());
                                        gr.upload_heights(&gpu.queue, &heights_init);
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
                                } else {
                                    // Unified progressive UI: use GPU raster for smooth rendering
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
                                                // Ensure GPU raster uses the same data as CPU overlays
                                                let verts = gpu_buf_mgr.get_depth_for_gpu(&world);
                                                rg.upload_inputs(&gpu.device, &gpu.queue, &u, face_ids.clone(), face_offs.clone(), &face_geom, verts, &face_edge_info, &face_perm_info);
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
                                                // println!("[viewer] raster(gpu) W={} H={} | F={} | verts={} | face_tbl={} | dispatch={}x{}", rw, rh, f_now, world.depth_m.len(), 20 * ((f_now + 1) * (f_now + 2) / 2), (rw + 7) / 8, (rh + 7) / 8);
                                            } else {
                                                // Always execute GPU raster when needed (data changed or texture missing)
                                                if ov.raster_dirty || ov.world_dirty || ov.color_dirty || ov.raster_tex_id.is_none() {
                                                    let mut dbg = if ov.gpu_dbg_wire { 1u32 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_face_tint { 1u32<<1 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_grid { 1u32<<2 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_tri_parity { 1u32<<3 } else { 0 };
                                                    dbg |= if ov.gpu_dbg_tri_index { 1u32<<5 } else { 0 };
                                                    dbg |= if ov.show_parity_heat || ov.export_parity_csv_requested { 1u32<<6 } else { 0 };
                                                    dbg |= if ov.force_cpu_face_pick { 1u32<<7 } else { 0 };
                                                    dbg |= if ov.dbg_cpu_bary_gpu_lattice { 1u32<<10 } else { 0 };
                                                    let u = raster_gpu::Uniforms { width: rw, height: rh, f: f_now, palette_mode: if ov.color_mode == 0 { 0 } else { 1 }, debug_flags: dbg, d_max: ov.hypso_d_max.max(1.0), h_max: ov.hypso_h_max.max(1.0), snowline: ov.hypso_snowline, eta_m: world.sea.eta_m, inv_dmax: 1.0f32/ov.hypso_d_max.max(1.0), inv_hmax: 1.0f32/ov.hypso_h_max.max(1.0) };
                                                    // Ensure GPU raster uses consistent data source
                                                    let verts = gpu_buf_mgr.get_depth_for_gpu(&world);
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
                                                    // Use the same vertex data for consistency
                                                    rg.write_vertex_values(&gpu.queue, &verts);
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
                                                                if false { if let Ok(mut fcsv) = std::fs::File::create("parity_roi.csv") {
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
                                                                } }
                                                                // Also dump only mismatching pixels across the full frame (capped)
                                                                if false { if let Ok(mut fmm) = std::fs::File::create("parity_mismatch.csv") {
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
                                                                } }
                                                            }
                                                        }
                                                    }
                                                    // Deferred CSV export (same pass) when requested
                                                    if ov.export_parity_csv_requested {
                                                        if let Some(buf) = rg.read_debug_face_tri(&gpu.device, &gpu.queue) {
                                                            // Export existing parity view for reference
                                                            let path = std::path::Path::new("parity_debug.csv");
                                                            if false { if let Ok(mut file) = std::fs::File::create(path) {
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
                                                            } }
                                                            println!("[parity] wrote parity_debug.csv ({}x{})", ov.raster_size.0, ov.raster_size.1);

                                                            // Export rollover probe CSV per checklist
                                                            let path_probe = std::path::Path::new("rollover_probe.csv");
                                                            if false { if let Ok(mut file2) = std::fs::File::create(path_probe) {
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
                                                            } }
                                                            println!("[probe] wrote rollover_probe.csv ({}x{})", ov.raster_size.0, ov.raster_size.1);
                                                        }
                                                        ov.export_parity_csv_requested = false;
                                                    }
                                                    if ov.raster_tex_id.is_none() { let tid = egui_renderer.register_native_texture(&gpu.device, &rg.out_view, wgpu::FilterMode::Linear); ov.raster_tex_id = Some(tid); }
                                                    ov.raster_dirty = false; ov.color_dirty = false; ov.last_raster_at = std::time::Instant::now();
                                                    // Don't clear world_dirty here - let overlays update first
                                                    // println!("[viewer] raster(gpu) W={} H={} | F={} | verts={} | face_tbl={} | dispatch={}x{}", rw, rh, f_now, world.depth_m.len(), 20 * ((f_now + 1) * (f_now + 2) / 2), (rw + 7) / 8, (rh + 7) / 8);
                                                }
                                            }
                                        }
                                        // Draw GPU raster and capture its actual rect for overlays
                                        let mut rect_img = rect;
                                        if let Some(tid) = ov.raster_tex_id {
                                            let avail = ui.available_size();
                                            let resp = ui.image(egui::load::SizedTexture::new(tid, avail));
                                            rect_img = resp.rect;
                                        }
                                        // Draw parity heat overlay as red dots if enabled
                                        if ov.show_parity_heat {
                                            if let Some(points) = &ov.parity_points {
                                                let rect_px = rect_img;
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
                                        // Ensure boundary classification reflects current velocities/plates
                                        if ov.world_dirty && (ov.show_bounds || ov.show_plate_type || ov.show_plate_adjacency || ov.show_triple_junctions) {
                                            world.boundaries = engine::boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
                                            ov.bounds_cache = None;
                                        }
                                        // Draw overlays on top (GPU raster provides the base elevation layer)
                                        overlay::draw_advanced_layers(ui, &painter, rect_img, &world, &world.grid, &mut ov);
                                        
                                        // Clear world_dirty flag after overlays have been updated
                                        ov.world_dirty = false;
                                    } else {
                                        // CPU fallback raster (only when GPU raster is disabled)
                                        // Keep boundary classification fresh when world changed
                                        if ov.world_dirty && (ov.show_bounds || ov.show_plate_type || ov.show_plate_adjacency || ov.show_triple_junctions) {
                                            world.boundaries = engine::boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
                                            ov.bounds_cache = None;
                                        }
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
                                        // Draw overlays on top (CPU raster provides the base elevation layer)
                                        overlay::draw_advanced_layers(ui, &painter, rect, &world, &world.grid, &mut ov);
                                        
                                        // Clear world_dirty flag after overlays have been updated
                                        ov.world_dirty = false;
                                    }
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
                                // If world changed, refresh vertex heights using elevation (η − depth)
                                if ov.world_dirty {
                                    let heights_now = gpu_buf_mgr.get_elevation_for_comparison(&world);
                                    gr.upload_heights(&gpu.queue, &heights_now);
                                }
                                gr.update_uniforms(
                                    &gpu.queue,
                                    view_proj,
                                    radius,
                                    exagger,
                                    dbg_flags,
                                    ov.hypso_d_max,
                                    ov.hypso_h_max,
                                    1.0f32 / get_phys_consts().r_earth_m as f32,
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
                        // Playback ticking - rate limited to prevent visual chaos
                        if ov.stepper.playing {
                            // Rate limit simulation steps to avoid overwhelming the viewer
                            let now = std::time::Instant::now();
                            let should_step = if let Some(last_step_time) = ov.last_sim_step_time {
                                now.duration_since(last_step_time).as_millis() >= 100 // Max 10 steps per second
                            } else {
                                true // First step
                            };
                            
                            if should_step && !sim_busy.load(Ordering::SeqCst) {
                                ov.last_sim_step_time = Some(now);
                                let cfg = engine::config::PipelineCfg {
                                    dt_myr: ov.sim_dt_myr.max(0.0),
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
                                let process_flags = ProcessFlags {
                                    enable_rigid_motion: ov.enable_rigid_motion,
                                    enable_subduction: ov.enable_subduction,
                                    enable_transforms: ov.enable_transforms,
                                    enable_flexure: ov.enable_flexure,
                                    enable_surface_processes: ov.enable_surface_processes,
                                    enable_isostasy: ov.enable_isostasy,
                                    enable_continental_buoyancy: ov.enable_continental_buoyancy,
                                    enable_orogeny: ov.enable_orogeny,
                                    enable_accretion: ov.enable_accretion,
                                    enable_rifting: ov.enable_rifting,
                                    enable_ridge_birth: ov.enable_ridge_birth,
                                    
                                    // Flexure backend configuration
                                    flexure_backend_cpu: ov.flexure_backend_cpu,
                                    flexure_gpu_levels: ov.flexure_gpu_levels,
                                    flexure_gpu_cycles: ov.flexure_gpu_cycles,
                                    
                                    // Unified cadence configuration
                                    cadence_config: {
                                        let mut config = engine::cadence_manager::CadenceConfig::new();
                                        config.set_cadence(engine::cadence_manager::ProcessType::RigidMotion, ov.cadence_rigid_motion);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Transforms, ov.cadence_transforms);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Subduction, ov.cadence_subduction);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Flexure, ov.cadence_flexure);
                                        config.set_cadence(engine::cadence_manager::ProcessType::SurfaceProcesses, ov.cadence_surface_processes);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Isostasy, ov.cadence_isostasy);
                                        config.set_cadence(engine::cadence_manager::ProcessType::ContinentalBuoyancy, ov.cadence_continental_buoyancy);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Orogeny, ov.cadence_orogeny);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Accretion, ov.cadence_accretion);
                                        config.set_cadence(engine::cadence_manager::ProcessType::Rifting, ov.cadence_rifting);
                                        config.set_cadence(engine::cadence_manager::ProcessType::RidgeBirth, ov.cadence_ridge_birth);
                                        config.set_cadence(engine::cadence_manager::ProcessType::ForceBalance, ov.cadence_force_balance);
                                        config
                                    },
                                };
                                let _ = tx_cmd.send(SimCommand::Step(cfg, process_flags));
                            }
                        }
                        // World updates now come through unified rx_world channel only
                        
                        // Atomic world updates - drain ALL pending updates before rendering (secondary location)
                        let mut latest_world_update: Option<WorldSnapshot> = None;
                        while let Ok(ws) = rx_world.try_recv() {
                            latest_world_update = Some(ws); // Keep only the most recent update
                        }
                        
                        // Apply the most recent world state atomically (if any)
                        if let Some(ws) = latest_world_update {
                            if world.depth_m.len() == ws.depth_m.len() { world.depth_m = ws.depth_m; }
                            if world.c.len() == ws.c.len() { world.c = ws.c; }
                            if world.th_c_m.len() == ws.th_c_m.len() { world.th_c_m = ws.th_c_m; }
                            world.sea.eta_m = ws.sea_eta_m;
                            if world.plates.plate_id.len() == ws.plate_id.len() { world.plates.plate_id = ws.plate_id; }
                            if world.plates.pole_axis.len() == ws.pole_axis.len() { world.plates.pole_axis = ws.pole_axis; }
                            if world.plates.omega_rad_yr.len() == ws.omega_rad_yr.len() { world.plates.omega_rad_yr = ws.omega_rad_yr; }
                            if world.v_en.len() == ws.v_en.len() { world.v_en = ws.v_en; }
                            if world.age_myr.len() == ws.age_myr.len() { world.age_myr = ws.age_myr; }
                            // Recompute boundaries with updated state
                            world.boundaries = engine::boundaries::Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, 0.005);
                            ov.world_dirty = true; ov.color_dirty = true; ov.raster_dirty = true;
                            // Invalidate overlay caches so they update with new world state
                            ov.plates_cache = None;
                            ov.vel_cache = None;
                            ov.bounds_cache = None;
                            ov.age_cache = None;
                            ov.bathy_cache = None;
                            ov.net_adj_cache = None;
                            ov.net_tj_cache = None;
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


fn elevation_curr_clone() -> Option<Vec<f32>> {
    // No longer using global elevation state - return None
    None
}
