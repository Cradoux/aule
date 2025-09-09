//! Flexible flexure backend management with automatic CPU/GPU fallback.
//!
//! This module provides a unified interface for flexure computation with
//! support for multiple backends and intelligent fallback strategies.

use crate::world::World;

/// Flexure computation backend selection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FlexureBackend {
    /// CPU-based Winkler foundation solver (always available)
    CpuWinkler,
    /// GPU-based multigrid solver with configurable parameters
    GpuMultigrid { 
        /// Number of multigrid levels
        levels: u32, 
        /// Number of V-cycles per solve
        cycles: u32 
    },
}

impl Default for FlexureBackend {
    fn default() -> Self {
        FlexureBackend::CpuWinkler
    }
}

/// Flexure computation manager with backend selection and fallback.
pub struct FlexureManager {
    /// Current backend configuration
    backend: FlexureBackend,
    /// Whether GPU backend is available
    gpu_available: bool,
    /// GPU backend parameters
    gpu_wj_omega: f32,
    /// Whether to subtract mean load for stability
    subtract_mean_load: bool,
}

impl FlexureManager {
    /// Create a new flexure manager with the specified backend.
    pub fn new(backend: FlexureBackend) -> Self {
        Self {
            backend,
            gpu_available: Self::check_gpu_availability(),
            gpu_wj_omega: 0.8, // Default weighted-Jacobi omega
            subtract_mean_load: true,
        }
    }

    /// Create a flexure manager that automatically selects the best available backend.
    pub fn auto() -> Self {
        let backend = if Self::check_gpu_availability() {
            FlexureBackend::GpuMultigrid { levels: 3, cycles: 2 }
        } else {
            FlexureBackend::CpuWinkler
        };
        Self::new(backend)
    }

    /// Get the current backend configuration.
    pub fn backend(&self) -> FlexureBackend {
        self.backend
    }

    /// Set the backend configuration.
    pub fn set_backend(&mut self, backend: FlexureBackend) {
        self.backend = backend;
    }

    /// Check if GPU backend is available.
    pub fn gpu_available(&self) -> bool {
        self.gpu_available
    }

    /// Set GPU parameters for multigrid solver.
    pub fn set_gpu_params(&mut self, wj_omega: f32, subtract_mean_load: bool) {
        self.gpu_wj_omega = wj_omega.clamp(0.4, 0.95);
        self.subtract_mean_load = subtract_mean_load;
    }

    /// Compute flexure using the configured backend with automatic fallback.
    /// 
    /// Returns the flexural deflection field and whether GPU was actually used.
    pub fn compute_flexure(
        &self,
        world: &World,
        f_load: &[f32],
    ) -> Result<(Vec<f32>, bool), FlexureError> {
        match self.backend {
            FlexureBackend::CpuWinkler => {
                let w_field = self.compute_cpu_winkler(world, f_load)?;
                Ok((w_field, false))
            }
            FlexureBackend::GpuMultigrid { levels, cycles } => {
                if self.gpu_available {
                    match self.compute_gpu_multigrid(world, f_load, levels, cycles) {
                        Ok(w_field) => Ok((w_field, true)),
                        Err(_) => {
                            // GPU failed, fallback to CPU
                            eprintln!("[flexure_manager] GPU flexure failed, falling back to CPU Winkler");
                            let w_field = self.compute_cpu_winkler(world, f_load)?;
                            Ok((w_field, false))
                        }
                    }
                } else {
                    // GPU not available, use CPU
                    let w_field = self.compute_cpu_winkler(world, f_load)?;
                    Ok((w_field, false))
                }
            }
        }
    }

    /// Check if GPU flexure is available on this system.
    fn check_gpu_availability() -> bool {
        // For now, always return true - actual availability will be tested at runtime
        // In the future, this could check for specific GPU features, memory, etc.
        true
    }

    /// CPU-based Winkler foundation flexure computation.
    fn compute_cpu_winkler(&self, world: &World, f_load: &[f32]) -> Result<Vec<f32>, FlexureError> {
        let n = world.grid.cells;
        if f_load.len() != n {
            return Err(FlexureError::InvalidLoadSize);
        }

        // Use the existing CPU flexure implementation
        let mut w_field = vec![0.0f32; n];
        
        // Apply mean load subtraction if enabled
        let load = if self.subtract_mean_load {
            let mean_load: f64 = f_load.iter().map(|&f| f as f64).sum::<f64>() / (n as f64);
            f_load.iter().map(|&f| f - mean_load as f32).collect()
        } else {
            f_load.to_vec()
        };

        // Compute flexure using existing Winkler solver
        // This is a simplified implementation - in practice, this would call
        // the existing flexure::compute_flexure_1d or similar
        for (i, &load_i) in load.iter().enumerate() {
            // Simplified elastic response (real implementation would use proper Green's functions)
            let te = world.te_m[i].max(1000.0); // Minimum 1km Te
            let area = world.area_m2[i];
            let response_factor = area / (te * te + 1e6); // Simplified elastic response
            w_field[i] = load_i * response_factor;
        }

        Ok(w_field)
    }

    /// GPU-based multigrid flexure computation.
    fn compute_gpu_multigrid(
        &self, 
        world: &World, 
        f_load: &[f32], 
        _levels: u32, 
        _cycles: u32
    ) -> Result<Vec<f32>, FlexureError> {
        let n = world.grid.cells;
        if f_load.len() != n {
            return Err(FlexureError::InvalidLoadSize);
        }

        // This would use the existing GPU flexure implementation from pipeline.rs
        // For now, we'll return an error to trigger CPU fallback
        // In a real implementation, this would call crate::flexure_gpu::FlexGpu

        // Apply mean load subtraction if enabled
        let _load = if self.subtract_mean_load {
            let mean_load: f64 = f_load.iter().map(|&f| f as f64).sum::<f64>() / (n as f64);
            f_load.iter().map(|&f| f - mean_load as f32).collect()
        } else {
            f_load.to_vec()
        };

        // Placeholder: Return error to trigger CPU fallback
        // Real implementation would use the existing GPU flexure code
        Err(FlexureError::GpuNotAvailable)
    }
}

/// Flexure computation errors.
#[derive(Debug, Clone)]
pub enum FlexureError {
    /// Invalid load field size
    InvalidLoadSize,
    /// GPU backend not available
    GpuNotAvailable,
    /// GPU computation failed
    GpuError(String),
    /// CPU computation failed
    CpuError(String),
}

impl std::fmt::Display for FlexureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FlexureError::InvalidLoadSize => write!(f, "Load field size does not match grid size"),
            FlexureError::GpuNotAvailable => write!(f, "GPU flexure backend not available"),
            FlexureError::GpuError(msg) => write!(f, "GPU flexure error: {}", msg),
            FlexureError::CpuError(msg) => write!(f, "CPU flexure error: {}", msg),
        }
    }
}

impl std::error::Error for FlexureError {}
