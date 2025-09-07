//! Configuration types for the unified pipeline.
//!
//! This module contains shared configuration structures used by both
//! the pipeline and world modules, breaking circular dependencies.

/// Pipeline execution mode that determines output behavior.
#[derive(Clone, Copy, Debug)]
pub enum PipelineMode {
    /// Realtime mode: solve eta separately, don't modify depth_m
    /// Renderer uses: elevation = eta - depth_m  
    Realtime { 
        /// If true, preserve depth_m buffer unchanged
        preserve_depth: bool 
    },
    /// Batch mode: write sea-level offset back into depth_m
    /// Renderer uses: elevation = -depth_m
    Batch { 
        /// If true, write sea-level offset back to depth_m
        write_back_sea_level: bool 
    },
}

impl Default for PipelineMode {
    fn default() -> Self {
        PipelineMode::Realtime { preserve_depth: true }
    }
}

/// Unified physics configuration that replaces both StepParams and PipelineCfg.
/// Simple and Advanced modes will use different presets of this configuration.
#[derive(Clone, Copy, Debug)]
pub struct PhysicsConfig {
    /// Time step in Myr
    pub dt_myr: f32,
    /// Number of steps to run per frame (viewer may loop externally as well)
    pub steps_per_frame: u32,
    
    // Process enables (independent flags)
    /// Enable rigid plate motion (advect plate_id, refresh velocities)
    pub enable_rigid_motion: bool,
    /// Enable subduction band edits
    pub enable_subduction: bool,
    /// Enable transform pull-apart/restraining bands  
    pub enable_transforms: bool,
    /// Enable elastic flexure response to current loads
    pub enable_flexure: bool,
    /// Enable surface processes (erosion, diffusion, sediment transport/deposition)
    pub enable_surface_processes: bool,
    /// Enable global sea level regulation to maintain reference ocean volume
    pub enable_isostasy: bool,
    /// Enable continental buoyancy response (separate from erosion)
    pub enable_continental_buoyancy: bool,
    /// Enable collision orogeny (C–C sutures)
    pub enable_orogeny: bool,
    /// Enable O–C accretion (arc/forearc growth)
    pub enable_accretion: bool,
    /// Enable continental rifting and passive margins
    pub enable_rifting: bool,
    /// Reset age along divergent boundaries (ridge births)
    pub enable_ridge_birth: bool,
    
    // Cadence controls (>=1); a stage runs when (step_idx+1) % cadence == 0
    /// Cadence for rigid motion updates
    pub cadence_rigid_motion: u32,
    /// Cadence for transform processing
    pub cadence_transforms: u32,
    /// Cadence for subduction processing
    pub cadence_subduction: u32,
    /// Cadence for flexure solving
    pub cadence_flexure: u32,
    /// Cadence for surface processes
    pub cadence_surface_processes: u32,
    /// Cadence for isostasy/sea-level regulation
    pub cadence_isostasy: u32,
    /// Cadence for continental buoyancy updates
    pub cadence_continental_buoyancy: u32,
    
    // Sea-level regulation
    /// Target land fraction (0..1) used by sea-level solve
    pub target_land_frac: f32,
    /// If true, keep eta fixed (skip sea-level regulation this tick)
    pub freeze_eta: bool,
    /// If true, auto re-baseline sea level after continents change
    pub auto_rebaseline_after_continents: bool,
    
    // Flexure backend selection
    /// Flexure computation backend
    pub flexure_backend: crate::flexure_manager::FlexureBackend,
    /// GPU flexure: weighted-Jacobi omega
    pub gpu_wj_omega: f32,
    /// Subtract mean load before solving (stability, removes rigid-body mode)
    pub subtract_mean_load: bool,
    
    // Legacy flexure configuration (deprecated, kept for compatibility)
    /// If true, attempt GPU flexure (experimental). Fallback to Winkler if unavailable
    pub use_gpu_flexure: bool,
    /// GPU flexure: number of multigrid levels
    pub gpu_flex_levels: u32,
    /// GPU flexure: V-cycles per cadence
    pub gpu_flex_cycles: u32,
    
    // Performance and stability
    /// Sub-steps for transforms per cadence (narrow operator stability)
    pub substeps_transforms: u32,
    /// Sub-steps for subduction band deltas per cadence
    pub substeps_subduction: u32,
    /// Surface processes sub-cycling count
    pub surf_subcycles: u32,
    
    // Logging and diagnostics
    /// If true, append a one-line mass budget log each step
    pub log_mass_budget: bool,
    
    // Surface processes parameters
    /// Parameters for surface processes (erosion, diffusion, sediment transport)
    pub surface_params: crate::surface::SurfaceParams,
    
    // Subduction parameters (viewer-controlled)
    /// Convergence threshold in m/yr used for band detection
    pub sub_tau_conv_m_per_yr: f32,
    /// Trench band half-width in km on subducting side
    pub sub_trench_half_width_km: f32,
    /// Arc band peak offset from trench hinge in km (overriding side)
    pub sub_arc_offset_km: f32,
    /// Arc band half-width in km
    pub sub_arc_half_width_km: f32,
    /// Back-arc band width in km (behind arc)
    pub sub_backarc_width_km: f32,
    /// Trench deepening magnitude in meters (positive deepens)
    pub sub_trench_deepen_m: f32,
    /// Arc uplift magnitude in meters (negative uplifts/shallows)
    pub sub_arc_uplift_m: f32,
    /// Back-arc uplift magnitude in meters (negative uplifts/shallows)
    pub sub_backarc_uplift_m: f32,
    /// Absolute rollback offset applied to overriding distances in meters
    pub sub_rollback_offset_m: f32,
    /// Rollback rate in km/Myr for time-progression (applied by viewer logic)
    pub sub_rollback_rate_km_per_myr: f32,
    /// If true, back-arc is deepened (extension mode) instead of uplifted
    pub sub_backarc_extension_mode: bool,
    /// Back-arc deepening magnitude in meters when in extension mode
    pub sub_backarc_extension_deepen_m: f32,
    /// Continental fraction threshold [0,1] for gating trench deepening
    pub sub_continent_c_min: f32,
}

impl PhysicsConfig {
    /// Configuration preset for Simple mode
    pub fn simple_mode() -> Self {
        Self {
            dt_myr: 1.0,
            steps_per_frame: 1,
            
            // Enable all processes every step
            enable_rigid_motion: true,
            enable_subduction: true,
            enable_transforms: true,
            enable_flexure: true,
            enable_surface_processes: false, // typically off in simple mode
            enable_isostasy: true,
            enable_continental_buoyancy: true,
            enable_orogeny: false,
            enable_accretion: false,
            enable_rifting: false,
            enable_ridge_birth: true,
            
            // Run every step (cadence = 1)
            cadence_rigid_motion: 1,
            cadence_transforms: 1,
            cadence_subduction: 1,
            cadence_flexure: 1,
            cadence_surface_processes: 1,
            cadence_isostasy: 1,
            cadence_continental_buoyancy: 1,
            
            // Sea-level
            target_land_frac: 0.3,
            freeze_eta: false,
            auto_rebaseline_after_continents: true,
            
            // Flexure backend settings
            flexure_backend: crate::flexure_manager::FlexureBackend::CpuWinkler,
            gpu_wj_omega: 0.8,
            subtract_mean_load: true,
            
            // Legacy GPU flexure settings (deprecated)
            use_gpu_flexure: false,
            gpu_flex_levels: 3,
            gpu_flex_cycles: 2,
            
            // Stability
            substeps_transforms: 4,
            substeps_subduction: 4,
            surf_subcycles: 1,
            
            // Logging
            log_mass_budget: false,
            
            // Default surface params (mostly unused in simple mode)
            surface_params: crate::surface::SurfaceParams::default(),
            
            // Default subduction parameters  
            sub_tau_conv_m_per_yr: 0.02,
            sub_trench_half_width_km: 50.0,
            sub_arc_offset_km: 150.0,
            sub_arc_half_width_km: 100.0,
            sub_backarc_width_km: 200.0,
            sub_trench_deepen_m: 1000.0,
            sub_arc_uplift_m: -500.0,
            sub_backarc_uplift_m: -200.0,
            sub_rollback_offset_m: 0.0,
            sub_rollback_rate_km_per_myr: 0.0,
            sub_backarc_extension_mode: false,
            sub_backarc_extension_deepen_m: 300.0,
            sub_continent_c_min: 0.5,
        }
    }
    
    /// Configuration preset for Advanced mode  
    pub fn advanced_mode() -> Self {
        Self {
            dt_myr: 0.5,
            steps_per_frame: 1,
            
            // Enable processes based on user toggles
            enable_rigid_motion: true,
            enable_subduction: true,
            enable_transforms: true,
            enable_flexure: true,
            enable_surface_processes: true,
            enable_isostasy: true,
            enable_continental_buoyancy: true,
            enable_orogeny: true,
            enable_accretion: true,
            enable_rifting: true,
            enable_ridge_birth: true,
            
            // Advanced mode uses different cadences for performance
            cadence_rigid_motion: 1,
            cadence_transforms: 2,
            cadence_subduction: 4,
            cadence_flexure: 1,
            cadence_surface_processes: 1,
            cadence_isostasy: 1,
            cadence_continental_buoyancy: 1,
            
            // Sea-level  
            target_land_frac: 0.3,
            freeze_eta: false,
            auto_rebaseline_after_continents: true,
            
            // Flexure backend settings (advanced mode can use GPU)
            flexure_backend: crate::flexure_manager::FlexureBackend::GpuMultigrid { levels: 3, cycles: 2 },
            gpu_wj_omega: 0.8,
            subtract_mean_load: true,
            
            // Legacy GPU flexure settings (deprecated)
            use_gpu_flexure: true, // Advanced mode tries GPU by default
            gpu_flex_levels: 3,
            gpu_flex_cycles: 2,
            
            // Stability
            substeps_transforms: 4,
            substeps_subduction: 4,
            surf_subcycles: 1,
            
            // Logging
            log_mass_budget: false,
            
            // Default surface params
            surface_params: crate::surface::SurfaceParams::default(),
            
            // Default subduction parameters
            sub_tau_conv_m_per_yr: 0.02,
            sub_trench_half_width_km: 50.0,
            sub_arc_offset_km: 150.0,
            sub_arc_half_width_km: 100.0,
            sub_backarc_width_km: 200.0,
            sub_trench_deepen_m: 1000.0,
            sub_arc_uplift_m: -500.0,
            sub_backarc_uplift_m: -200.0,
            sub_rollback_offset_m: 0.0,
            sub_rollback_rate_km_per_myr: 0.0,
            sub_backarc_extension_mode: false,
            sub_backarc_extension_deepen_m: 300.0,
            sub_continent_c_min: 0.5,
        }
    }
    
    /// Convert to legacy PipelineCfg format for compatibility
    pub fn to_pipeline_cfg(&self) -> PipelineCfg {
        PipelineCfg {
            dt_myr: self.dt_myr,
            steps_per_frame: self.steps_per_frame,
            enable_flexure: self.enable_flexure,
            enable_erosion: self.enable_surface_processes,
            target_land_frac: self.target_land_frac,
            freeze_eta: self.freeze_eta,
            log_mass_budget: self.log_mass_budget,
            enable_subduction: self.enable_subduction,
            enable_rigid_motion: self.enable_rigid_motion,
            cadence_trf_every: self.cadence_transforms,
            cadence_sub_every: self.cadence_subduction,
            cadence_flx_every: self.cadence_flexure,
            cadence_sea_every: self.cadence_isostasy,
            cadence_surf_every: self.cadence_surface_processes,
            substeps_transforms: self.substeps_transforms,
            substeps_subduction: self.substeps_subduction,
            use_gpu_flexure: self.use_gpu_flexure,
            gpu_flex_levels: self.gpu_flex_levels,
            gpu_flex_cycles: self.gpu_flex_cycles,
            gpu_wj_omega: self.gpu_wj_omega,
            subtract_mean_load: self.subtract_mean_load,
            sub_tau_conv_m_per_yr: self.sub_tau_conv_m_per_yr,
            sub_trench_half_width_km: self.sub_trench_half_width_km,
            sub_arc_offset_km: self.sub_arc_offset_km,
            sub_arc_half_width_km: self.sub_arc_half_width_km,
            sub_backarc_width_km: self.sub_backarc_width_km,
            sub_trench_deepen_m: self.sub_trench_deepen_m,
            sub_arc_uplift_m: self.sub_arc_uplift_m,
            sub_backarc_uplift_m: self.sub_backarc_uplift_m,
            sub_rollback_offset_m: self.sub_rollback_offset_m,
            sub_rollback_rate_km_per_myr: self.sub_rollback_rate_km_per_myr,
            sub_backarc_extension_mode: self.sub_backarc_extension_mode,
            sub_backarc_extension_deepen_m: self.sub_backarc_extension_deepen_m,
            sub_continent_c_min: self.sub_continent_c_min,
            surf_k_stream: self.surface_params.k_stream,
            surf_m_exp: self.surface_params.m_exp,
            surf_n_exp: self.surface_params.n_exp,
            surf_k_diff: self.surface_params.k_diff,
            surf_k_tr: self.surface_params.k_tr,
            surf_p_exp: self.surface_params.p_exp,
            surf_q_exp: self.surface_params.q_exp,
            surf_rho_sed: self.surface_params.rho_sed,
            surf_min_slope: self.surface_params.min_slope,
            surf_subcycles: self.surf_subcycles,
            surf_couple_flexure: self.surface_params.couple_flexure,
            cadence_spawn_plate_every: 0, // disabled by default
            cadence_retire_plate_every: 0, // disabled by default  
            cadence_force_balance_every: 0, // disabled by default
            fb_gain: 1.0e-12,
            fb_damp_per_myr: 0.2,
            fb_k_conv: 1.0,
            fb_k_div: 0.5,
            fb_k_trans: 0.1,
            fb_max_domega: 5.0e-9,
            fb_max_omega: 2.0e-7,
        }
    }
    
    /// Convert to legacy StepParams format for compatibility
    pub fn to_step_params(&self) -> StepParams {
        StepParams {
            dt_myr: self.dt_myr as f64,
            do_flexure: self.enable_flexure,
            do_isostasy: self.enable_isostasy,
            do_transforms: self.enable_transforms,
            do_subduction: self.enable_subduction,
            do_continents: self.enable_continental_buoyancy,
            do_ridge_birth: self.enable_ridge_birth,
            auto_rebaseline_after_continents: self.auto_rebaseline_after_continents,
            do_rigid_motion: self.enable_rigid_motion,
            do_orogeny: self.enable_orogeny,
            do_accretion: self.enable_accretion,
            do_rifting: self.enable_rifting,
            do_surface: self.enable_surface_processes,
            surface_params: self.surface_params,
            advection_every: self.cadence_rigid_motion,
            transforms_every: self.cadence_transforms,
            subduction_every: self.cadence_subduction,
            flexure_every: self.cadence_flexure,
            sea_every: self.cadence_isostasy,
            do_advection: self.enable_rigid_motion,
            do_sea: self.enable_isostasy,
        }
    }
}

/// Parameters that control a single evolution step (legacy world::step_once).
#[derive(Clone, Copy, Debug)]
pub struct StepParams {
    /// Time step in Myr
    pub dt_myr: f64,
    /// Apply elastic flexure response to current loads
    pub do_flexure: bool,
    /// Adjust global sea level to maintain reference ocean volume
    pub do_isostasy: bool,
    /// Apply transform pull-apart/restraining bands
    pub do_transforms: bool,
    /// Apply subduction trench/arc/backarc edits
    pub do_subduction: bool,
    /// Advect C/th_c and apply continental uplift to depth
    pub do_continents: bool,
    /// Reset age along divergent boundaries (ridge births)
    pub do_ridge_birth: bool,
    /// If true, auto re-baseline sea level after continents change.
    pub auto_rebaseline_after_continents: bool,
    /// Enable rigid plate motion (advect plate_id and update velocities)
    pub do_rigid_motion: bool,
    /// Enable collision orogeny (C–C sutures)
    pub do_orogeny: bool,
    /// Enable O–C accretion (arc/forearc growth)
    pub do_accretion: bool,
    /// Enable continental rifting and passive margins
    pub do_rifting: bool,
    /// Enable surface processes (erosion, diffusion, sediment transport/deposition)
    pub do_surface: bool,
    /// Parameter set for surface processes.
    pub surface_params: crate::surface::SurfaceParams,
    /// Cadence: run advection every N steps (>=1). When 1, runs each step.
    pub advection_every: u32,
    /// Cadence for transforms
    pub transforms_every: u32,
    /// Cadence for subduction
    pub subduction_every: u32,
    /// Cadence for flexure
    pub flexure_every: u32,
    /// Cadence for sea-level/isostasy
    pub sea_every: u32,
    /// Gate advection explicitly (combined with cadence)
    pub do_advection: bool,
    /// Gate sea-level explicitly (combined with cadence)
    pub do_sea: bool,
}

/// Pipeline configuration used by Simple and Advanced modes (from pipeline.rs).
#[derive(Clone, Copy, Debug)]
pub struct PipelineCfg {
    /// Time step in Myr.
    pub dt_myr: f32,
    /// Number of steps to run per frame (viewer may loop externally as well).
    pub steps_per_frame: u32,
    /// Enable flexure solve stage.
    pub enable_flexure: bool,
    /// Enable erosion/diffusion stage.
    pub enable_erosion: bool,
    /// Target land fraction (0..1) used by sea-level solve.
    pub target_land_frac: f32,
    /// If true, keep eta fixed (skip sea-level regulation this tick).
    pub freeze_eta: bool,
    /// If true, append a one-line mass budget log each step.
    pub log_mass_budget: bool,
    /// Enable subduction band edits (disable to isolate physics noise)
    pub enable_subduction: bool,
    /// Enable rigid plate motion (advect `plate_id`, refresh velocities)
    pub enable_rigid_motion: bool,
    /// Cadence controls (>=1); a stage runs when `(step_idx+1) % cadence == 0`.
    /// Transformer bands cadence in steps.
    pub cadence_trf_every: u32,
    /// Subduction bands cadence in steps.
    pub cadence_sub_every: u32,
    /// Flexure cadence in steps.
    pub cadence_flx_every: u32,
    /// Sea-level (eta) solve cadence in steps.
    pub cadence_sea_every: u32,
    /// Surface processes cadence in steps.
    pub cadence_surf_every: u32,
    /// Sub-steps for transforms per cadence (narrow operator stability)
    pub substeps_transforms: u32,
    /// Sub-steps for subduction band deltas per cadence
    pub substeps_subduction: u32,
    /// If true, attempt GPU flexure (experimental). Fallback to Winkler if unavailable.
    pub use_gpu_flexure: bool,
    /// GPU flexure: number of multigrid levels.
    pub gpu_flex_levels: u32,
    /// GPU flexure: V-cycles per cadence.
    pub gpu_flex_cycles: u32,
    /// GPU flexure: weighted-Jacobi omega.
    pub gpu_wj_omega: f32,
    /// Subtract mean load before solving (stability, removes rigid-body mode).
    pub subtract_mean_load: bool,
    // Subduction tunables (viewer-controlled)
    /// Convergence threshold in m/yr used for band detection.
    pub sub_tau_conv_m_per_yr: f32,
    /// Trench band half-width in km on subducting side.
    pub sub_trench_half_width_km: f32,
    /// Arc band peak offset from trench hinge in km (overriding side).
    pub sub_arc_offset_km: f32,
    /// Arc band half-width in km.
    pub sub_arc_half_width_km: f32,
    /// Back-arc band width in km (behind arc).
    pub sub_backarc_width_km: f32,
    /// Trench deepening magnitude in meters (positive deepens).
    pub sub_trench_deepen_m: f32,
    /// Arc uplift magnitude in meters (negative uplifts/shallows).
    pub sub_arc_uplift_m: f32,
    /// Back-arc uplift magnitude in meters (negative uplifts/shallows).
    pub sub_backarc_uplift_m: f32,
    /// Absolute rollback offset applied to overriding distances in meters.
    pub sub_rollback_offset_m: f32,
    /// Rollback rate in km/Myr for time-progression (applied by viewer logic).
    pub sub_rollback_rate_km_per_myr: f32,
    /// If true, back-arc is deepened (extension mode) instead of uplifted.
    pub sub_backarc_extension_mode: bool,
    /// Back-arc deepening magnitude in meters when in extension mode.
    pub sub_backarc_extension_deepen_m: f32,
    /// Continental fraction threshold [0,1] for gating trench deepening.
    pub sub_continent_c_min: f32,

    // Surface processes parameters (viewer-controlled)
    /// Stream-power coefficient for fluvial incision (yr⁻¹ m^(1-m) s^m units folded)
    pub surf_k_stream: f32,
    /// Stream-power m exponent
    pub surf_m_exp: f32,
    /// Stream-power n exponent
    pub surf_n_exp: f32,
    /// Hillslope diffusion coefficient κ (m²/yr)
    pub surf_k_diff: f32,
    /// Sediment transport coefficient (nondimensional scaling)
    pub surf_k_tr: f32,
    /// Transport-limited p exponent
    pub surf_p_exp: f32,
    /// Transport-limited q exponent
    pub surf_q_exp: f32,
    /// Sediment density (kg/m³)
    pub surf_rho_sed: f32,
    /// Minimum slope for slope-dependent processes (dimensionless)
    pub surf_min_slope: f32,
    /// Surface processes sub-cycling count
    pub surf_subcycles: u32,
    /// If true, couple surface processes with flexural response
    pub surf_couple_flexure: bool,

    // Advanced cadence controls
    /// Spawn plate cadence (0 = disabled)
    pub cadence_spawn_plate_every: u32,
    /// Retire plate cadence (0 = disabled)
    pub cadence_retire_plate_every: u32,
    /// Force balance update cadence (0 = disabled)
    pub cadence_force_balance_every: u32,
    /// Force balance gain parameter
    pub fb_gain: f32,
    /// Force balance damping per Myr
    pub fb_damp_per_myr: f32,
    /// Force balance convergent coefficient
    pub fb_k_conv: f32,
    /// Force balance divergent coefficient
    pub fb_k_div: f32,
    /// Force balance transform coefficient
    pub fb_k_trans: f32,
    /// Max absolute change in omega per step
    pub fb_max_domega: f32,
    /// Max absolute |ω| clamp.
    pub fb_max_omega: f32,
}

impl PhysicsConfig {
    /// Convert from legacy PipelineCfg for backward compatibility during migration.
    pub fn from_pipeline_cfg(cfg: PipelineCfg) -> Self {
        Self {
            dt_myr: cfg.dt_myr,
            enable_rigid_motion: cfg.enable_rigid_motion,
            enable_flexure: cfg.enable_flexure,
            enable_surface_processes: cfg.enable_erosion,
            enable_continental_buoyancy: true, // Always enabled in unified pipeline
            enable_subduction: cfg.enable_subduction,
            target_land_frac: cfg.target_land_frac,
            freeze_eta: cfg.freeze_eta,
            steps_per_frame: cfg.steps_per_frame,
            log_mass_budget: cfg.log_mass_budget,
            
            // Flexure settings
            use_gpu_flexure: cfg.use_gpu_flexure,
            gpu_flex_levels: cfg.gpu_flex_levels,
            gpu_flex_cycles: cfg.gpu_flex_cycles,
            gpu_wj_omega: cfg.gpu_wj_omega,
            subtract_mean_load: cfg.subtract_mean_load,
            
            // Stability
            substeps_transforms: cfg.substeps_transforms,
            substeps_subduction: cfg.substeps_subduction,
            surf_subcycles: cfg.surf_subcycles,
            
            // Copy surface params
            surface_params: crate::surface::SurfaceParams {
                k_stream: cfg.surf_k_stream,
                m_exp: cfg.surf_m_exp,
                n_exp: cfg.surf_n_exp,
                k_diff: cfg.surf_k_diff,
                k_tr: cfg.surf_k_tr,
                p_exp: cfg.surf_p_exp,
                q_exp: cfg.surf_q_exp,
                rho_sed: cfg.surf_rho_sed,
                min_slope: cfg.surf_min_slope,
                subcycles: cfg.surf_subcycles,
                couple_flexure: cfg.surf_couple_flexure,
            },
            
            // Copy subduction params  
            sub_tau_conv_m_per_yr: cfg.sub_tau_conv_m_per_yr,
            sub_trench_half_width_km: cfg.sub_trench_half_width_km,
            sub_arc_offset_km: cfg.sub_arc_offset_km,
            sub_arc_half_width_km: cfg.sub_arc_half_width_km,
            sub_backarc_width_km: cfg.sub_backarc_width_km,
            sub_trench_deepen_m: cfg.sub_trench_deepen_m,
            sub_arc_uplift_m: cfg.sub_arc_uplift_m,
            sub_backarc_uplift_m: cfg.sub_backarc_uplift_m,
            sub_rollback_offset_m: cfg.sub_rollback_offset_m,
            sub_rollback_rate_km_per_myr: cfg.sub_rollback_rate_km_per_myr,
            sub_backarc_extension_mode: cfg.sub_backarc_extension_mode,
            sub_backarc_extension_deepen_m: cfg.sub_backarc_extension_deepen_m,
            sub_continent_c_min: cfg.sub_continent_c_min,
            
            // Force balance - not directly available in PipelineCfg, use defaults
            ..Self::simple_mode()
        }
    }
}
