//! Direct GPU buffer access with no caching or conversion layers.
//!
//! This module provides direct access to world data for GPU operations,
//! eliminating redundant caches and conversion steps.

use engine::world::World;

/// Direct GPU buffer manager - no caching, no conversion.
pub struct GpuBufferManager;

impl GpuBufferManager {
    /// Create a new GPU buffer manager.
    pub fn new() -> Self {
        Self
    }

    /// Get depth data for GPU upload.
    /// 
    /// Returns direct reference to world.depth_m (positive down) for GPU shaders.
    /// The GPU shader will compute elevation as eta_m - depth.
    pub fn get_depth_for_gpu<'a>(&self, world: &'a World) -> &'a [f32] {
        &world.depth_m
    }

    /// Get elevation data for CPU/GPU comparison.
    /// 
    /// Returns elevation values (positive up) computed directly from world state.
    pub fn get_elevation_for_comparison(&self, world: &World) -> Vec<f32> {
        world.depth_m.iter().map(|&d| world.sea.eta_m - d).collect()
    }

    /// No-op for compatibility - no cache to invalidate.
    pub fn invalidate_cache(&mut self) {
        // No cache to invalidate
    }
}

impl Default for GpuBufferManager {
    fn default() -> Self {
        Self::new()
    }
}
