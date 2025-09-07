//! GPU buffer management with consistent elevation data synchronization.
//!
//! This module provides a centralized manager for GPU buffer operations
//! to ensure consistency between CPU and GPU elevation data.

use crate::get_elevation_state;
use engine::world::World;

/// Centralized GPU buffer manager for consistent elevation data.
pub struct GpuBufferManager {
    /// Cached elevation data to avoid redundant computations
    cached_elevation: Option<Vec<f32>>,
    /// Whether the cached elevation is valid
    cache_valid: bool,
}

impl GpuBufferManager {
    /// Create a new GPU buffer manager.
    pub fn new() -> Self {
        Self {
            cached_elevation: None,
            cache_valid: false,
        }
    }

    /// Get consistent elevation data for GPU upload.
    /// 
    /// This method ensures all GPU operations use the same elevation calculation:
    /// - Priority 1: Use simulation thread elevation from ElevationState
    /// - Priority 2: Fallback to world.depth_m converted to elevation
    /// 
    /// Returns depth values (positive down) for GPU shaders.
    pub fn get_depth_for_gpu(&mut self, world: &World) -> Vec<f32> {
        if !self.cache_valid {
            let elevation = if let Some(elev) = get_elevation_state().get_clone() {
                // Use elevation from simulation thread (most up-to-date)
                elev.iter().map(|&z| world.sea.eta_m - z).collect()
            } else {
                // Fallback to world depth_m (static world state)
                world.depth_m.clone()
            };
            
            self.cached_elevation = Some(elevation);
            self.cache_valid = true;
        }
        
        self.cached_elevation.clone().unwrap_or_else(|| world.depth_m.clone())
    }

    /// Get consistent elevation data for CPU/GPU comparison.
    /// 
    /// Returns elevation values (positive up) for consistency checks.
    pub fn get_elevation_for_comparison(&mut self, world: &World) -> Vec<f32> {
        if let Some(elev) = get_elevation_state().get_clone() {
            elev
        } else {
            // Convert depth to elevation: elevation = eta - depth
            world.depth_m.iter().map(|&d| world.sea.eta_m - d).collect()
        }
    }

    /// Invalidate the cache when world state changes.
    pub fn invalidate_cache(&mut self) {
        self.cache_valid = false;
        self.cached_elevation = None;
    }

    /// Check if the cache is valid.
    pub fn is_cache_valid(&self) -> bool {
        self.cache_valid
    }
}

impl Default for GpuBufferManager {
    fn default() -> Self {
        Self::new()
    }
}
