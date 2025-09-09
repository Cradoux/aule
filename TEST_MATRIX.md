# Aulé Test Matrix

## Quick Acceptance Tests

### Continental Motion Alignment
**Objective**: Verify continental centroids move with plate motion
**Script**: 
```rust
#[test]
fn test_continental_centroid_alignment() {
    let mut world = setup_test_world_with_continents();
    let initial_centroids = compute_continental_centroids(&world);
    
    // Run 10 Myr with known plate motion
    let dt_myr = 1.0;
    for _ in 0..10 {
        unified_pipeline.step(&mut world, PipelineMode::Batch);
    }
    
    let final_centroids = compute_continental_centroids(&world);
    let predicted_centroids = predict_centroids_from_plates(&initial_centroids, &world.plates, 10.0);
    
    for (final, predicted) in final_centroids.iter().zip(predicted_centroids.iter()) {
        let error_deg = haversine_distance_deg(*final, *predicted);
        assert!(error_deg < 0.25, "Continental centroid error: {:.3}°", error_deg);
    }
}
```
**Acceptance**: <0.25° error for all continental centroids
**Frequency**: Every PR that touches advection

### Buoyancy Independence Test
**Objective**: Verify buoyancy operates without erosion enabled
**Script**:
```rust
#[test]
fn test_buoyancy_without_erosion() {
    let mut world = setup_world_with_thick_continents();
    let config = PhysicsConfig {
        enable_surface_processes: false,  // erosion OFF
        enable_continental_buoyancy: true,  // buoyancy ON
        ..Default::default()
    };
    
    let initial_depth = world.depth_m.clone();
    unified_pipeline.step(&mut world, config, PipelineMode::Batch);
    
    let max_depth_change = world.depth_m.iter()
        .zip(initial_depth.iter())
        .map(|(new, old)| (new - old).abs())
        .fold(0.0f32, f32::max);
    
    assert!(max_depth_change > 1.0, "Expected >1m depth change from buoyancy, got {:.3}m", max_depth_change);
}
```
**Acceptance**: `dZ_step_max > 1 m` when buoyancy enabled, regardless of erosion flag
**Frequency**: Every PR that touches process flags

### Renderer Buffer Consistency
**Objective**: Verify renderer uses same elevation buffer in both modes
**Script**:
```rust
#[test]
fn test_renderer_buffer_consistency() {
    let mut world = setup_identical_test_world();
    
    // Simple mode
    let config_simple = PhysicsConfig::simple_mode();
    let mut pipeline_simple = UnifiedPipeline::new(config_simple);
    pipeline_simple.step(&mut world, PipelineMode::Realtime);
    let elevation_simple = pipeline_simple.elevation_manager().elevation().to_vec();
    
    // Reset world
    world = setup_identical_test_world();
    
    // Advanced mode  
    let config_advanced = PhysicsConfig::advanced_mode();
    let mut pipeline_advanced = UnifiedPipeline::new(config_advanced);
    pipeline_advanced.step(&mut world, PipelineMode::Batch);
    let elevation_advanced = pipeline_advanced.elevation_manager().elevation().to_vec();
    
    // Compare elevations
    let max_diff = elevation_simple.iter()
        .zip(elevation_advanced.iter())
        .map(|(s, a)| (s - a).abs())
        .fold(0.0f32, f32::max);
    
    assert!(max_diff < 0.1, "Mode elevation difference: {:.3}m", max_diff);
}
```
**Acceptance**: Elevation buffers identical between modes within 0.1m
**Frequency**: Every PR that touches pipeline or elevation management

### Sea-Level Freeze Test
**Objective**: Verify sea-level freeze affects only flooding, not elevation
**Script**:
```rust
#[test]
fn test_sea_level_freeze_isolation() {
    let mut world = setup_test_world();
    let config = PhysicsConfig {
        freeze_sea_level: false,
        enable_subduction: true,
        ..Default::default()
    };
    
    // Step with sea-level regulation
    let mut pipeline = UnifiedPipeline::new(config);
    pipeline.step(&mut world, PipelineMode::Realtime);
    let depth_unfrozen = world.depth_m.clone();
    let eta_unfrozen = world.sea.eta_m;
    
    // Reset and step with frozen sea-level
    world = setup_test_world();
    config.freeze_sea_level = true;
    pipeline = UnifiedPipeline::new(config);
    pipeline.step(&mut world, PipelineMode::Realtime);
    let depth_frozen = world.depth_m.clone();
    let eta_frozen = world.sea.eta_m;
    
    // Depth should be identical (tectonic processes unaffected)
    let max_depth_diff = depth_unfrozen.iter()
        .zip(depth_frozen.iter())
        .map(|(u, f)| (u - f).abs())
        .fold(0.0f32, f32::max);
    
    assert!(max_depth_diff < 0.01, "Depth changed with sea-level freeze: {:.3}m", max_depth_diff);
    
    // Eta should be different (regulation vs freeze)
    let eta_diff = (eta_unfrozen - eta_frozen).abs();
    assert!(eta_diff > 0.1, "Sea-level should differ when frozen: {:.3}m", eta_diff);
}
```
**Acceptance**: Depth values unchanged, only eta differs when sea-level frozen
**Frequency**: Every PR that touches sea-level or isostasy

## Performance Sanity Tests

### Frame Time Benchmark
**Objective**: Ensure frame times remain acceptable
**Script**:
```rust
#[test]
fn benchmark_frame_time() {
    let mut world = setup_large_world(frequency = 256); // ~65k cells
    let config = PhysicsConfig::simple_mode();
    let mut pipeline = UnifiedPipeline::new(config);
    
    let iterations = 100;
    let start = std::time::Instant::now();
    
    for _ in 0..iterations {
        pipeline.step(&mut world, PipelineMode::Realtime);
    }
    
    let elapsed = start.elapsed();
    let avg_frame_ms = elapsed.as_secs_f64() * 1000.0 / (iterations as f64);
    
    println!("Average frame time: {:.1}ms", avg_frame_ms);
    assert!(avg_frame_ms < 16.7, "Frame time too slow: {:.1}ms", avg_frame_ms); // 60 FPS target
}
```
**Acceptance**: <16.7ms average frame time at F=256 on target hardware
**Frequency**: Every PR that touches core pipeline

### Memory Usage Test
**Objective**: Detect memory leaks and excessive allocation
**Script**:
```rust
#[test]
fn test_memory_stability() {
    let mut world = setup_test_world();
    let config = PhysicsConfig::default();
    let mut pipeline = UnifiedPipeline::new(config);
    
    let initial_memory = get_process_memory_mb();
    
    // Run many steps
    for _ in 0..1000 {
        pipeline.step(&mut world, PipelineMode::Realtime);
    }
    
    let final_memory = get_process_memory_mb();
    let memory_growth = final_memory - initial_memory;
    
    assert!(memory_growth < 50.0, "Memory leak detected: +{:.1}MB", memory_growth);
}
```
**Acceptance**: <50MB memory growth over 1000 steps
**Frequency**: Weekly regression test

### GPU/CPU Transfer Efficiency
**Objective**: Minimize unnecessary GPU↔CPU transfers
**Script**:
```rust
#[test]
fn test_gpu_transfer_efficiency() {
    let mut world = setup_test_world();
    let config = PhysicsConfig {
        flexure_backend: FlexureBackend::GpuMultigrid { levels: 3, cycles: 2 },
        ..Default::default()
    };
    let mut pipeline = UnifiedPipeline::new(config);
    
    let transfer_counter = GpuTransferCounter::new();
    
    for _ in 0..10 {
        pipeline.step(&mut world, PipelineMode::Realtime);
    }
    
    let transfers = transfer_counter.get_stats();
    println!("GPU transfers: {} up, {} down, {:.1}MB total", 
             transfers.uploads, transfers.downloads, transfers.total_mb);
    
    // Should not download elevation data every frame
    assert!(transfers.downloads < 5, "Too many GPU downloads: {}", transfers.downloads);
}
```
**Acceptance**: <5 GPU downloads per 10 frames
**Frequency**: Every PR that touches GPU code

## Integration Test Suite

### End-to-End Simulation Test
**Objective**: Run complete 100 Myr simulation without crashes
**Script**:
```rust
#[test]
#[ignore] // Long-running test
fn test_long_simulation() {
    let mut world = setup_earth_like_world();
    let config = PhysicsConfig::advanced_mode();
    let mut pipeline = UnifiedPipeline::new(config);
    
    let target_time_myr = 100.0;
    let dt_myr = 0.5;
    let steps = (target_time_myr / dt_myr) as usize;
    
    for i in 0..steps {
        pipeline.step(&mut world, PipelineMode::Batch);
        
        // Sanity checks every 10 Myr
        if i % 20 == 0 {
            validate_world_sanity(&world);
        }
    }
    
    assert_eq!(world.clock.t_myr, target_time_myr);
    validate_final_state(&world);
}
```
**Acceptance**: Completes without panics, final state passes validation
**Frequency**: Nightly regression test

### Mode Equivalence Test
**Objective**: Verify Simple and Advanced modes produce equivalent results
**Script**:
```rust
#[test]
fn test_mode_equivalence() {
    let world_template = setup_deterministic_world(seed = 12345);
    
    // Run Simple mode
    let mut world_simple = world_template.clone();
    let config_simple = PhysicsConfig::simple_mode();
    let mut pipeline_simple = UnifiedPipeline::new(config_simple);
    
    for _ in 0..20 { // 20 Myr
        pipeline_simple.step(&mut world_simple, PipelineMode::Realtime);
    }
    
    // Run Advanced mode
    let mut world_advanced = world_template.clone();
    let config_advanced = PhysicsConfig::advanced_mode();
    let mut pipeline_advanced = UnifiedPipeline::new(config_advanced);
    
    for _ in 0..20 {
        pipeline_advanced.step(&mut world_advanced, PipelineMode::Batch);
    }
    
    // Compare key metrics
    let simple_stats = compute_world_stats(&world_simple);
    let advanced_stats = compute_world_stats(&world_advanced);
    
    assert_relative_eq!(simple_stats.land_fraction, advanced_stats.land_fraction, epsilon = 0.02);
    assert_relative_eq!(simple_stats.mean_elevation, advanced_stats.mean_elevation, epsilon = 10.0);
    assert_relative_eq!(simple_stats.continental_area, advanced_stats.continental_area, epsilon = 0.05);
}
```
**Acceptance**: Key statistics within 2-5% between modes
**Frequency**: Every PR that changes pipeline behavior

## Regression Test Matrix

### Visual Regression Tests
**Objective**: Detect visual changes in rendered output
**Tools**: Image comparison with tolerance
**Frequency**: Every PR
**Acceptance**: <1% pixel difference in hypsometric rendering

### Determinism Tests
**Objective**: Same inputs produce identical outputs
**Script**: Run identical simulation twice, compare all outputs bit-for-bit
**Frequency**: Every PR that touches core simulation
**Acceptance**: Bit-perfect reproduction

### Boundary Condition Tests
**Objective**: Handle edge cases gracefully
**Cases**:
- Zero-sized continents
- Single-cell plates  
- Extreme time steps (0.001 Myr, 10 Myr)
- Grid frequencies (F=16 to F=1024)
- All processes disabled
**Acceptance**: No panics, reasonable fallback behavior

## Automated Test Infrastructure

### Continuous Integration Pipeline
```yaml
name: Aulë Test Suite
on: [push, pull_request]
jobs:
  quick-tests:
    runs-on: [ubuntu-latest, windows-latest, macos-latest]
    steps:
      - run: cargo test --lib
      - run: cargo test test_continental_centroid_alignment
      - run: cargo test test_buoyancy_without_erosion  
      - run: cargo test test_renderer_buffer_consistency
      - run: cargo test benchmark_frame_time
      
  integration-tests:
    runs-on: ubuntu-latest
    steps:
      - run: cargo test --test integration -- --ignored
      - run: cargo test test_mode_equivalence
      
  nightly-tests:
    runs-on: ubuntu-latest
    if: github.event_name == 'schedule'
    steps:
      - run: cargo test test_long_simulation -- --ignored
      - run: cargo test test_memory_stability
```

### Test Data Management
- **Golden datasets**: Stored snapshots of known-good simulation states
- **Regression baselines**: Performance and memory benchmarks
- **Visual references**: PNG outputs for image comparison
- **Determinism seeds**: Fixed random seeds for reproducible tests

### Test Reporting
- **Coverage reports**: Code coverage for all modules
- **Performance trends**: Frame time and memory usage over time  
- **Visual diffs**: Highlighted changes in rendered output
- **Failure analysis**: Automatic bisection for regression identification

## Test Environment Requirements

### Hardware Specifications
- **Minimum**: Intel i5-8400 / AMD Ryzen 5 2600, 16GB RAM, GTX 1060 / RX 580
- **Target**: Intel i7-13700K / AMD Ryzen 7 7700X, 32GB RAM, RTX 4080 / RX 7800XT
- **CI**: GitHub Actions standard runners (2-core, 7GB RAM)

### Software Dependencies
- **Rust**: Latest stable (1.75+)
- **GPU Drivers**: Latest stable for target hardware
- **Test Tools**: `cargo-nextest`, `cargo-criterion`, custom visual diff tools

### Test Data Size Limits
- **Unit tests**: <1MB total test data
- **Integration tests**: <100MB including golden datasets
- **Long-running tests**: <1GB temporary data, cleaned after run

This test matrix ensures comprehensive validation of the unification plan while maintaining reasonable CI performance and test maintainability.
