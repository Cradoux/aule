//! Unified pipeline that replaces both world::step_once and pipeline::step_full.
//!
//! This module provides a single execution path that both Simple and Advanced modes
//! use, with mode differences expressed as configuration rather than separate code paths.

use crate::{
    config::{PhysicsConfig, PipelineMode}, 
    elevation_manager::ElevationManager, 
    pipeline::SurfaceFields, world::World
};

/// Result of a unified pipeline step.
#[derive(Clone, Copy, Debug)]
pub struct StepResult {
    /// Simulation time after the step (Myr)
    pub t_myr: f64,
    /// Time step size used (Myr)
    pub dt_myr: f64,
    /// Number of divergent boundary edges
    pub div_count: u32,
    /// Number of convergent boundary edges  
    pub conv_count: u32,
    /// Number of transform boundary edges
    pub trans_count: u32,
    /// Area-weighted mean of continental fraction C (0..1)
    pub c_bar: f64,
    /// Flexure residual if applied; otherwise negative
    pub flex_residual: f32,
    /// Sea level eta after regulation
    pub eta_m: f32,
    /// Performance timings (ms)
    pub timings: StepTimings,
}

/// Performance timing breakdown for a pipeline step.
#[derive(Clone, Copy, Debug, Default)]
pub struct StepTimings {
    /// Time spent on boundary classification (ms)
    pub boundaries_ms: f64,
    /// Time spent on kinematics/advection (ms)
    pub kinematics_ms: f64,
    /// Time spent on ridge processing (ms)
    pub ridge_ms: f64,
    /// Time spent on subduction processing (ms)
    pub subduction_ms: f64,
    /// Time spent on transform processing (ms)
    pub transforms_ms: f64,
    /// Time spent on continental buoyancy (ms)
    pub buoyancy_ms: f64,
    /// Time spent on flexure solving (ms)
    pub flexure_ms: f64,
    /// Time spent on surface processes (ms)
    pub surface_ms: f64,
    /// Time spent on sea-level regulation (ms)
    pub sea_level_ms: f64,
    /// Total time for entire step (ms)
    pub total_ms: f64,
}

/// Unified pipeline executor that handles both Simple and Advanced mode behavior.
pub struct UnifiedPipeline {
    config: PhysicsConfig,
    /// Centralized elevation management with intelligent caching
    elevation_mgr: Option<ElevationManager>,
}

impl UnifiedPipeline {
    /// Create a new unified pipeline with the given configuration.
    pub fn new(config: PhysicsConfig) -> Self {
        Self {
            config,
            elevation_mgr: None,
        }
    }
    
    /// Update the pipeline configuration.
    pub fn update_config(&mut self, config: PhysicsConfig) {
        self.config = config;
        self.invalidate_elevation_cache();
    }
    
    /// Get current configuration.
    pub fn config(&self) -> &PhysicsConfig {
        &self.config
    }
    
    /// Execute one physics step with configurable output mode.
    pub fn step(&mut self, world: &mut World, mode: PipelineMode) -> StepResult {
        let _start_time = std::time::Instant::now();
        let timings = StepTimings::default();
        
        // Invalidate elevation cache since world will change
        self.invalidate_elevation_cache();
        
        // UNIFIED APPROACH: Always use the same pipeline logic, but handle sea-level output differently
        let mut elevation_temp = vec![0.0f32; world.grid.cells];
        let mut eta_temp = world.sea.eta_m;
        let surf = SurfaceFields {
            elevation_m: &mut elevation_temp,
            eta_m: &mut eta_temp,
        };
        
        // Execute the unified pipeline step directly (without deprecated step_full)
        self.execute_physics_step(world, surf);
        
        // Handle sea-level output based on mode
        match mode {
            PipelineMode::Realtime { preserve_depth: _ } => {
                // Realtime mode: Keep eta separate, don't modify depth_m
                world.sea.eta_m = eta_temp;
                // Elevation manager will handle caching automatically
            }
            
            PipelineMode::Batch { write_back_sea_level } => {
                // Batch mode: Write sea-level offset back into depth_m if requested
                if write_back_sea_level {
                    let eta_delta = eta_temp - world.sea.eta_m;
                    for depth in world.depth_m.iter_mut() {
                        *depth += eta_delta;
                    }
                    // Reset eta to 0 since offset is now in depth_m
                    world.sea.eta_m = 0.0;
                } else {
                    // Just update eta without writing back
                    world.sea.eta_m = eta_temp;
                }
                // Cache elevation for renderer (always eta - depth_m for consistency)
                // Elevation manager will compute elevation on demand
            }
        }
        
        // Extract results from world state
        StepResult {
            t_myr: world.clock.t_myr,
            dt_myr: self.config.dt_myr as f64,
            div_count: world.boundaries.stats.divergent,
            conv_count: world.boundaries.stats.convergent,
            trans_count: world.boundaries.stats.transform,
            c_bar: {
                let mut c_sum = 0.0f64;
                let mut area_sum = 0.0f64;
                for (c, area) in world.c.iter().zip(world.area_m2.iter()) {
                    c_sum += (*c as f64) * (*area as f64);
                    area_sum += *area as f64;
                }
                if area_sum > 0.0 { c_sum / area_sum } else { 0.0 }
            },
            flex_residual: world.last_flex_residual,
            eta_m: world.sea.eta_m,
            timings,
        }
    }
    
    /// Get current elevation field (positive up).
    /// Returns cached elevation if valid, otherwise computes from depth_m and eta_m.
    pub fn elevation(&mut self, world: &World) -> &[f32] {
        // Initialize elevation manager if needed
        if self.elevation_mgr.is_none() {
            self.elevation_mgr = Some(ElevationManager::new(
                world.depth_m.clone(),
                world.sea.eta_m
            ));
        }
        
        // Update elevation manager with current world state and get elevation
        if let Some(ref mut mgr) = self.elevation_mgr {
            mgr.update_depth_and_eta(world.depth_m.clone(), world.sea.eta_m);
            mgr.get_elevation()
        } else {
            unreachable!("Elevation manager should be initialized above")
        }
    }
    
    /// Check if elevation data has changed since last access.
    pub fn is_elevation_dirty(&self) -> bool {
        self.elevation_mgr.as_ref()
            .map(|mgr| !mgr.is_cache_valid())
            .unwrap_or(true)
    }
    
    /// Mark elevation as needing recomputation.
    pub fn invalidate_elevation_cache(&mut self) {
        if let Some(ref mut mgr) = self.elevation_mgr {
            mgr.invalidate();
        }
    }
    
    /// Get a cloned copy of the elevation field.
    pub fn elevation_clone(&mut self, world: &World) -> Vec<f32> {
        self.elevation(world).to_vec()
    }
    
    /// Execute the core physics step with unified configuration.
    /// This replaces the deprecated pipeline::step_full function.
    fn execute_physics_step(&mut self, world: &mut World, surf: SurfaceFields) {
        // Convert PhysicsConfig to PipelineCfg for the execution
        let pipeline_cfg = self.config.to_pipeline_cfg();
        
        // Call the core physics implementation (will be fully internalized in future cleanup)
        crate::pipeline::step_full_internal(world, surf, pipeline_cfg);
    }
}

impl Default for UnifiedPipeline {
    fn default() -> Self {
        Self::new(PhysicsConfig::simple_mode())
    }
}
