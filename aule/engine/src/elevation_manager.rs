//! Centralized elevation buffer management.
//!
//! This module provides a single source of truth for elevation data with
//! intelligent caching and invalidation to eliminate redundant computations.

/// Centralized elevation buffer manager with caching and invalidation.
pub struct ElevationManager {
    /// Primary depth field (positive down, meters)
    depth_m: Vec<f32>,
    /// Sea level offset (meters above geoid)
    eta_m: f32,
    /// Cached elevation field (positive up, meters) 
    cached_elevation: Option<Vec<f32>>,
    /// Whether the cached elevation is valid
    cache_valid: bool,
}

impl ElevationManager {
    /// Create a new elevation manager with initial depth field.
    pub fn new(depth_m: Vec<f32>, eta_m: f32) -> Self {
        Self {
            depth_m,
            eta_m,
            cached_elevation: None,
            cache_valid: false,
        }
    }

    /// Get the current elevation field (positive up), computing if needed.
    /// This is the primary interface for accessing elevation data.
    pub fn get_elevation(&mut self) -> &[f32] {
        if !self.cache_valid || self.cached_elevation.is_none() {
            self.recompute_elevation();
        }
        self.cached_elevation.as_ref().unwrap()
    }

    /// Get a cloned copy of the elevation field, computing if needed.
    pub fn get_elevation_clone(&mut self) -> Vec<f32> {
        self.get_elevation().to_vec()
    }

    /// Update the depth field and invalidate cache.
    pub fn update_depth(&mut self, new_depth: Vec<f32>) {
        if self.depth_m.len() != new_depth.len() {
            // Resize cached elevation if grid size changed
            self.cached_elevation = None;
        }
        self.depth_m = new_depth;
        self.invalidate();
    }

    /// Update the sea level eta and invalidate cache.
    pub fn update_eta(&mut self, new_eta: f32) {
        if (self.eta_m - new_eta).abs() > 1e-6 {
            self.eta_m = new_eta;
            self.invalidate();
        }
    }

    /// Update both depth and eta simultaneously (more efficient).
    pub fn update_depth_and_eta(&mut self, new_depth: Vec<f32>, new_eta: f32) {
        if self.depth_m.len() != new_depth.len() {
            self.cached_elevation = None;
        }
        self.depth_m = new_depth;
        self.eta_m = new_eta;
        self.invalidate();
    }

    /// Mark the elevation cache as invalid, forcing recomputation on next access.
    pub fn invalidate(&mut self) {
        self.cache_valid = false;
    }

    /// Check if the elevation cache is currently valid.
    pub fn is_cache_valid(&self) -> bool {
        self.cache_valid && self.cached_elevation.is_some()
    }

    /// Get direct access to depth field (for compatibility).
    pub fn depth_m(&self) -> &[f32] {
        &self.depth_m
    }

    /// Get current eta value.
    pub fn eta_m(&self) -> f32 {
        self.eta_m
    }

    /// Get the number of cells in the elevation field.
    pub fn len(&self) -> usize {
        self.depth_m.len()
    }

    /// Check if the elevation manager is empty.
    pub fn is_empty(&self) -> bool {
        self.depth_m.is_empty()
    }

    /// Recompute the elevation field from depth and eta.
    fn recompute_elevation(&mut self) {
        let n = self.depth_m.len();
        
        // Initialize or resize cached elevation if needed
        match &mut self.cached_elevation {
            Some(elevation) if elevation.len() == n => {
                // Reuse existing buffer
                for (i, &depth) in self.depth_m.iter().enumerate() {
                    elevation[i] = self.eta_m - depth;
                }
            }
            _ => {
                // Create new buffer
                self.cached_elevation = Some(
                    self.depth_m.iter()
                        .map(|&depth| self.eta_m - depth)
                        .collect()
                );
            }
        }
        
        self.cache_valid = true;
    }
}

/// Trait for objects that can provide elevation data.
pub trait ElevationProvider {
    /// Get current elevation field (positive up)
    fn elevation(&mut self) -> &[f32];
    
    /// Check if elevation data has changed since last access
    fn is_dirty(&self) -> bool;
    
    /// Mark elevation as needing recomputation
    fn invalidate(&mut self);
}

impl ElevationProvider for ElevationManager {
    fn elevation(&mut self) -> &[f32] {
        self.get_elevation()
    }
    
    fn is_dirty(&self) -> bool {
        !self.is_cache_valid()
    }
    
    fn invalidate(&mut self) {
        self.invalidate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_elevation_computation() {
        let depth = vec![100.0, -50.0, 200.0]; // Mixed depths
        let eta = 10.0;
        let mut mgr = ElevationManager::new(depth, eta);
        
        let elevation = mgr.get_elevation();
        assert_eq!(elevation[0], -90.0);  // eta - depth = 10 - 100 = -90
        assert_eq!(elevation[1], 60.0);   // eta - depth = 10 - (-50) = 60  
        assert_eq!(elevation[2], -190.0); // eta - depth = 10 - 200 = -190
    }

    #[test]
    fn test_cache_invalidation() {
        let depth = vec![100.0, 200.0];
        let eta = 0.0;
        let mut mgr = ElevationManager::new(depth, eta);
        
        // First access computes cache
        let _elev1 = mgr.get_elevation();
        assert!(mgr.is_cache_valid());
        
        // Update depth invalidates cache
        mgr.update_depth(vec![150.0, 250.0]);
        assert!(!mgr.is_cache_valid());
        
        // Access recomputes cache
        {
            let elev2 = mgr.get_elevation();
            assert_eq!(elev2[0], -150.0); // New depth reflected
        }
        assert!(mgr.is_cache_valid()); // Check cache validity after using the elevation
    }

    #[test]
    fn test_eta_update() {
        let depth = vec![100.0];
        let mut mgr = ElevationManager::new(depth, 0.0);
        
        let elev1 = mgr.get_elevation();
        assert_eq!(elev1[0], -100.0);
        
        mgr.update_eta(50.0);
        let elev2 = mgr.get_elevation();
        assert_eq!(elev2[0], -50.0); // eta change reflected
    }
}
