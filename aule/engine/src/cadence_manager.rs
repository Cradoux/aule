//! Unified cadence management system for all physics processes.
//!
//! This module provides a centralized way to configure and manage
//! the execution cadence (frequency) of different physics processes.

use std::collections::HashMap;

/// Physics process identifiers for cadence management.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ProcessType {
    /// Rigid plate motion and kinematics
    RigidMotion,
    /// Transform fault processing
    Transforms,
    /// Subduction band processing  
    Subduction,
    /// Flexural isostasy computation
    Flexure,
    /// Surface processes (erosion, diffusion)
    SurfaceProcesses,
    /// Isostatic sea-level regulation
    Isostasy,
    /// Continental buoyancy updates
    ContinentalBuoyancy,
    /// Orogenic processes
    Orogeny,
    /// Accretion processes
    Accretion,
    /// Rifting processes
    Rifting,
    /// Ridge birth processes
    RidgeBirth,
    /// Force balance updates
    ForceBalance,
    /// Plate spawning
    PlateSpawning,
    /// Plate retirement
    PlateRetirement,
}

impl ProcessType {
    /// Get all available process types.
    pub fn all() -> Vec<ProcessType> {
        vec![
            ProcessType::RigidMotion,
            ProcessType::Transforms,
            ProcessType::Subduction,
            ProcessType::Flexure,
            ProcessType::SurfaceProcesses,
            ProcessType::Isostasy,
            ProcessType::ContinentalBuoyancy,
            ProcessType::Orogeny,
            ProcessType::Accretion,
            ProcessType::Rifting,
            ProcessType::RidgeBirth,
            ProcessType::ForceBalance,
            ProcessType::PlateSpawning,
            ProcessType::PlateRetirement,
        ]
    }

    /// Get a human-readable name for the process.
    pub fn display_name(&self) -> &'static str {
        match self {
            ProcessType::RigidMotion => "Rigid Motion",
            ProcessType::Transforms => "Transforms",
            ProcessType::Subduction => "Subduction",
            ProcessType::Flexure => "Flexure",
            ProcessType::SurfaceProcesses => "Surface Processes",
            ProcessType::Isostasy => "Isostasy",
            ProcessType::ContinentalBuoyancy => "Continental Buoyancy",
            ProcessType::Orogeny => "Orogeny",
            ProcessType::Accretion => "Accretion",
            ProcessType::Rifting => "Rifting",
            ProcessType::RidgeBirth => "Ridge Birth",
            ProcessType::ForceBalance => "Force Balance",
            ProcessType::PlateSpawning => "Plate Spawning",
            ProcessType::PlateRetirement => "Plate Retirement",
        }
    }

    /// Get default cadence value for this process.
    pub fn default_cadence(&self) -> u32 {
        match self {
            // Core processes run every step in simple mode
            ProcessType::RigidMotion => 1,
            ProcessType::Flexure => 1,
            ProcessType::Isostasy => 1,
            ProcessType::ContinentalBuoyancy => 1,
            ProcessType::SurfaceProcesses => 1,
            
            // More expensive processes run less frequently
            ProcessType::Transforms => 2,
            ProcessType::Subduction => 4,
            ProcessType::Orogeny => 4,
            ProcessType::Accretion => 8,
            ProcessType::Rifting => 8,
            ProcessType::RidgeBirth => 8,
            
            // Infrequent processes
            ProcessType::ForceBalance => 8,
            ProcessType::PlateSpawning => 0, // Disabled by default
            ProcessType::PlateRetirement => 0, // Disabled by default
        }
    }

    /// Get advanced mode cadence value for this process.
    pub fn advanced_cadence(&self) -> u32 {
        match self {
            // Core processes still run frequently in advanced mode
            ProcessType::RigidMotion => 1,
            ProcessType::Flexure => 1,
            ProcessType::Isostasy => 1,
            ProcessType::ContinentalBuoyancy => 1,
            
            // Advanced users may want different performance trade-offs
            ProcessType::Transforms => 2,
            ProcessType::Subduction => 4,
            ProcessType::SurfaceProcesses => 2,
            ProcessType::Orogeny => 4,
            ProcessType::Accretion => 8,
            ProcessType::Rifting => 8,
            ProcessType::RidgeBirth => 8,
            
            // Infrequent processes
            ProcessType::ForceBalance => 8,
            ProcessType::PlateSpawning => 0, // Disabled by default
            ProcessType::PlateRetirement => 0, // Disabled by default
        }
    }

    /// Check if this process is considered "core" (should run frequently).
    pub fn is_core_process(&self) -> bool {
        matches!(
            self,
            ProcessType::RigidMotion
                | ProcessType::Flexure
                | ProcessType::Isostasy
                | ProcessType::ContinentalBuoyancy
        )
    }

    /// Check if this process is computationally expensive.
    pub fn is_expensive(&self) -> bool {
        matches!(
            self,
            ProcessType::Subduction
                | ProcessType::SurfaceProcesses
                | ProcessType::Orogeny
                | ProcessType::Accretion
                | ProcessType::Rifting
        )
    }
}

/// Cadence configuration for physics processes.
#[derive(Debug, Clone)]
pub struct CadenceConfig {
    /// Per-process cadence values (steps between executions)
    cadences: HashMap<ProcessType, u32>,
}

impl CadenceConfig {
    /// Create a new cadence configuration with default values.
    pub fn new() -> Self {
        let mut cadences = HashMap::new();
        for process_type in ProcessType::all() {
            cadences.insert(process_type, process_type.default_cadence());
        }
        Self { cadences }
    }

    /// Create a cadence configuration optimized for simple mode.
    pub fn simple_mode() -> Self {
        let mut cadences = HashMap::new();
        for process_type in ProcessType::all() {
            cadences.insert(process_type, process_type.default_cadence());
        }
        Self { cadences }
    }

    /// Create a cadence configuration optimized for advanced mode.
    pub fn advanced_mode() -> Self {
        let mut cadences = HashMap::new();
        for process_type in ProcessType::all() {
            cadences.insert(process_type, process_type.advanced_cadence());
        }
        Self { cadences }
    }

    /// Get cadence for a specific process type.
    pub fn get_cadence(&self, process_type: ProcessType) -> u32 {
        self.cadences
            .get(&process_type)
            .copied()
            .unwrap_or_else(|| process_type.default_cadence())
    }

    /// Set cadence for a specific process type.
    pub fn set_cadence(&mut self, process_type: ProcessType, cadence: u32) {
        self.cadences.insert(process_type, cadence);
    }

    /// Check if a process should execute on the given step.
    pub fn should_execute(&self, process_type: ProcessType, step: u64) -> bool {
        let cadence = self.get_cadence(process_type);
        if cadence == 0 {
            false // Process disabled
        } else {
            step % (cadence as u64) == 0
        }
    }

    /// Apply performance preset: "Balanced", "Performance", "Quality".
    pub fn apply_preset(&mut self, preset: &str) {
        match preset {
            "Performance" => {
                // Reduce cadences for better performance
                self.set_cadence(ProcessType::Transforms, 4);
                self.set_cadence(ProcessType::Subduction, 8);
                self.set_cadence(ProcessType::SurfaceProcesses, 4);
                self.set_cadence(ProcessType::Orogeny, 8);
                self.set_cadence(ProcessType::Accretion, 16);
                self.set_cadence(ProcessType::Rifting, 16);
            }
            "Quality" => {
                // Increase cadences for better quality
                self.set_cadence(ProcessType::Transforms, 1);
                self.set_cadence(ProcessType::Subduction, 2);
                self.set_cadence(ProcessType::SurfaceProcesses, 1);
                self.set_cadence(ProcessType::Orogeny, 2);
                self.set_cadence(ProcessType::Accretion, 4);
                self.set_cadence(ProcessType::Rifting, 4);
            }
            "Balanced" | _ => {
                // Use default/advanced cadences
                for process_type in ProcessType::all() {
                    self.set_cadence(process_type, process_type.advanced_cadence());
                }
            }
        }
    }

    /// Get all cadences as a vector for UI display.
    pub fn get_all_cadences(&self) -> Vec<(ProcessType, u32)> {
        let mut result = Vec::new();
        for process_type in ProcessType::all() {
            result.push((process_type, self.get_cadence(process_type)));
        }
        result
    }

    /// Update from legacy PipelineCfg for compatibility.
    pub fn update_from_pipeline_cfg(&mut self, cfg: &crate::config::PipelineCfg) {
        self.set_cadence(ProcessType::Transforms, cfg.cadence_trf_every);
        self.set_cadence(ProcessType::Subduction, cfg.cadence_sub_every);
        self.set_cadence(ProcessType::Flexure, cfg.cadence_flx_every);
        self.set_cadence(ProcessType::Isostasy, cfg.cadence_sea_every);
        self.set_cadence(ProcessType::SurfaceProcesses, cfg.cadence_surf_every);
    }
}

impl Default for CadenceConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// Cadence manager that tracks execution state across steps.
pub struct CadenceManager {
    /// Current cadence configuration
    config: CadenceConfig,
    /// Current step number
    current_step: u64,
}

impl CadenceManager {
    /// Create a new cadence manager with the given configuration.
    pub fn new(config: CadenceConfig) -> Self {
        Self {
            config,
            current_step: 0,
        }
    }

    /// Get the current cadence configuration.
    pub fn config(&self) -> &CadenceConfig {
        &self.config
    }

    /// Get a mutable reference to the cadence configuration.
    pub fn config_mut(&mut self) -> &mut CadenceConfig {
        &mut self.config
    }

    /// Set a new cadence configuration.
    pub fn set_config(&mut self, config: CadenceConfig) {
        self.config = config;
    }

    /// Advance to the next step.
    pub fn advance_step(&mut self) {
        self.current_step += 1;
    }

    /// Get the current step number.
    pub fn current_step(&self) -> u64 {
        self.current_step
    }

    /// Check if a process should execute on the current step.
    pub fn should_execute(&self, process_type: ProcessType) -> bool {
        self.config.should_execute(process_type, self.current_step)
    }

    /// Reset step counter (for simulation restart).
    pub fn reset(&mut self) {
        self.current_step = 0;
    }
}
