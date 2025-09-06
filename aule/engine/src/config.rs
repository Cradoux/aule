//! Configuration types for the unified pipeline.
//!
//! This module contains shared configuration structures used by both
//! the pipeline and world modules, breaking circular dependencies.

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
