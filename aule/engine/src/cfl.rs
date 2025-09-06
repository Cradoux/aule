//! CFL (Courant-Friedrichs-Lewy) limiter for numerical stability in tectonic processes.
//!
//! This module provides a unified approach to limiting displacement rates based on
//! characteristic widths and time steps to maintain numerical stability.
//!
//! For any process that moves material a distance Δ = u * dt across a characteristic width W,
//! compute C = |Δ| / W and scale the update by min(1, CFL_MAX / C) with CFL_MAX = 0.3.

/// CFL limiter configuration.
#[derive(Clone, Copy, Debug)]
pub struct CflConfig {
    /// Maximum allowed CFL number (default: 0.3).
    pub max_cfl: f64,
    /// Whether to enable debug logging of CFL violations.
    pub debug_log: bool,
}

impl Default for CflConfig {
    fn default() -> Self {
        Self { max_cfl: 0.3, debug_log: false }
    }
}

/// CFL limiter result containing the scaled displacement and diagnostic information.
#[derive(Clone, Copy, Debug)]
pub struct CflResult {
    /// Scaled displacement in meters.
    pub scaled_displacement_m: f64,
    /// Raw CFL number before scaling.
    pub raw_cfl: f64,
    /// CFL number after scaling (raw_cfl * scale_factor).
    pub capped_cfl: f64,
    /// Applied scaling factor (≤ 1.0).
    pub scale_factor: f64,
    /// Whether the displacement was scaled due to CFL violation.
    pub was_scaled: bool,
}

/// Limit displacement based on CFL condition.
///
/// # Arguments
/// * `raw_displacement_m` - Raw displacement in meters for this time step
/// * `width_m` - Characteristic width of the process in meters
/// * `dt_myr` - Time step in millions of years
/// * `config` - CFL configuration
///
/// # Returns
/// Scaled displacement that satisfies CFL ≤ max_cfl.
///
/// # Formula
/// CFL = |displacement| / width
/// scale_factor = min(1.0, max_cfl / CFL)
/// scaled_displacement = raw_displacement * scale_factor
pub fn limit(raw_displacement_m: f64, width_m: f64, dt_myr: f64, config: CflConfig) -> CflResult {
    // Convert dt from Myr to years for completeness (not needed when passing per-step displacement)
    let _dt_yr = dt_myr * 1.0e6;

    // Compute raw CFL: displacement / width
    let raw_cfl = if width_m > 0.0 { raw_displacement_m.abs() / width_m } else { f64::INFINITY };

    // Apply CFL limiting
    let scale_factor = if raw_cfl > config.max_cfl { config.max_cfl / raw_cfl } else { 1.0 };

    let scaled_displacement_m = raw_displacement_m * scale_factor;
    let was_scaled = scale_factor < 1.0;

    // Debug logging if enabled
    if config.debug_log && was_scaled {
        println!(
            "[CFL] limiting displacement: raw={:.3} m, width={:.0} m, CFL={:.3} → {:.3} (scale={:.3})",
            raw_displacement_m, width_m, raw_cfl, raw_cfl * scale_factor, scale_factor
        );
    }

    CflResult {
        scaled_displacement_m,
        raw_cfl,
        capped_cfl: raw_cfl * scale_factor,
        scale_factor,
        was_scaled,
    }
}

/// Convenience wrapper matching the RL-1 spec signature.
///
/// Returns only the scaled displacement in meters.
pub fn limit_simple(raw_displacement_m: f64, width_m: f64, dt_myr: f64, cfl_max: f64) -> f64 {
    let res = limit(
        raw_displacement_m,
        width_m,
        dt_myr,
        CflConfig { max_cfl: cfl_max, debug_log: false },
    );
    res.scaled_displacement_m
}

/// CFL statistics for a process over multiple cells.
#[derive(Debug)]
pub struct CflStats {
    /// Number of cells processed.
    pub cells_processed: u32,
    /// Number of cells that required CFL limiting.
    pub cells_scaled: u32,
    /// Maximum capped CFL encountered (after scaling).
    pub max_cfl: f64,
    /// Mean capped CFL across all cells.
    pub mean_cfl: f64,
    /// Minimum scale factor applied.
    pub min_scale_factor: f64,
}

impl CflStats {
    /// Update statistics with a new CFL result.
    pub fn update(&mut self, result: &CflResult) {
        self.cells_processed += 1;
        if result.was_scaled {
            self.cells_scaled += 1;
        }
        // Track capped CFL (after scaling) so acceptance (max ≤ CFL_MAX + 5%) is visible in logs
        self.max_cfl = self.max_cfl.max(result.capped_cfl);
        self.min_scale_factor = self.min_scale_factor.min(result.scale_factor);
        // Update mean capped CFL (running average)
        let n = self.cells_processed as f64;
        self.mean_cfl = (self.mean_cfl * (n - 1.0) + result.capped_cfl) / n;
    }

    /// Print summary statistics.
    pub fn print_summary(&self, process_name: &str) {
        if self.cells_processed > 0 {
            let scaled_pct = (self.cells_scaled as f64 / self.cells_processed as f64) * 100.0;
            println!(
                "[CFL] {}: {} cells, {:.1}% scaled, CFL max/mean={:.3}/{:.3}, min_scale={:.3}",
                process_name,
                self.cells_processed,
                scaled_pct,
                self.max_cfl,
                self.mean_cfl,
                self.min_scale_factor
            );
        }
    }
}

impl Default for CflStats {
    fn default() -> Self {
        Self {
            cells_processed: 0,
            cells_scaled: 0,
            max_cfl: 0.0,
            mean_cfl: 0.0,
            min_scale_factor: 1.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cfl_limit_below_threshold() {
        let config = CflConfig::default();
        let result = limit(100.0, 1000.0, 0.001, config);

        assert!(!result.was_scaled);
        assert_eq!(result.scale_factor, 1.0);
        assert_eq!(result.scaled_displacement_m, 100.0);
        assert_eq!(result.raw_cfl, 0.1);
    }

    #[test]
    fn test_cfl_limit_above_threshold() {
        let config = CflConfig::default();
        let result = limit(500.0, 1000.0, 0.001, config);

        assert!(result.was_scaled);
        assert_eq!(result.scale_factor, 0.3 / 0.5);
        assert_eq!(result.scaled_displacement_m, 300.0);
        assert_eq!(result.raw_cfl, 0.5);
    }

    #[test]
    fn test_cfl_limit_zero_width() {
        let config = CflConfig::default();
        let result = limit(100.0, 0.0, 0.001, config);

        assert!(result.was_scaled);
        assert_eq!(result.scale_factor, 0.0);
        assert_eq!(result.scaled_displacement_m, 0.0);
        assert!(result.raw_cfl.is_infinite());
    }

    #[test]
    fn test_cfl_stats() {
        let mut stats = CflStats::default();

        let result1 = limit(100.0, 1000.0, 0.001, CflConfig::default());
        let result2 = limit(500.0, 1000.0, 0.001, CflConfig::default());

        stats.update(&result1);
        stats.update(&result2);

        assert_eq!(stats.cells_processed, 2);
        assert_eq!(stats.cells_scaled, 1);
        // With capped CFL tracked, max_cfl reflects scaled value ≤ max_cfl
        assert!(stats.max_cfl <= CflConfig::default().max_cfl + 1e-12);
        assert!((stats.mean_cfl - 0.2).abs() < 1e-12 || (stats.mean_cfl - 0.3).abs() < 1e-12);
        assert_eq!(stats.min_scale_factor, 0.6);
    }

    #[test]
    fn test_cfl_property_bounded_for_various_inputs() {
        let cfg = CflConfig::default();
        let widths = [0.0, 1.0, 10.0, 1_000.0, 50_000.0];
        let raws = [-10_000.0, -1_000.0, -10.0, -1.0, 0.0, 1.0, 10.0, 1_000.0, 10_000.0];
        let dts = [0.0, 1e-6, 0.1, 1.0];
        for &w in &widths {
            for &raw in &raws {
                for &dt in &dts {
                    let res = limit(raw, w, dt, cfg);
                    if w > 0.0 {
                        let c = (res.scaled_displacement_m.abs()) / w;
                        assert!(
                            c <= cfg.max_cfl + 1e-12,
                            "CFL bound violated: c={} > max={} (raw={}, w={}, dt={})",
                            c,
                            cfg.max_cfl,
                            raw,
                            w,
                            dt
                        );
                    } else {
                        // w == 0 handled as infinite raw_cfl and zeroed displacement
                        assert_eq!(res.scaled_displacement_m, 0.0);
                        assert!(res.raw_cfl.is_infinite());
                    }
                }
            }
        }
    }

    #[test]
    fn test_cfl_monotonicity_abs_scaled_vs_abs_raw() {
        let cfg = CflConfig::default();
        let w = 1_000.0; // meters
        let dt = 0.5; // Myr (not used in formula, retained for signature)
                      // Monotonicity in |raw|: test on non-negative raw so |raw| increases with raw
        let inputs = [0.0, 0.01, 0.1, 0.5, 1.0, 10.0, 100.0, 500.0, 1_000.0, 10_000.0, 1_000_000.0];
        let mut prev_abs_scaled = 0.0;
        for &raw in &inputs {
            let res = limit(raw, w, dt, cfg);
            let abs_scaled = res.scaled_displacement_m.abs();
            // Non-decreasing with respect to |raw|
            assert!(
                abs_scaled + 1e-12 >= prev_abs_scaled,
                "monotonicity violated: prev={} curr={} (raw={})",
                prev_abs_scaled,
                abs_scaled,
                raw
            );
            prev_abs_scaled = abs_scaled;
            // Upper bound plateau ≈ cfg.max_cfl * w
            assert!(
                abs_scaled <= cfg.max_cfl * w + 1e-9,
                "scaled exceeds plateau: {} > {}",
                abs_scaled,
                cfg.max_cfl * w
            );
        }
        // Check strong-violation case hits the plateau within tolerance
        let res_big = limit(1.0e12, w, dt, cfg);
        let plateau = cfg.max_cfl * w;
        assert!((res_big.scaled_displacement_m.abs() - plateau).abs() <= 1e-9);
    }
}
