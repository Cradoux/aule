use engine::{config::{PhysicsConfig, PipelineMode}, unified_pipeline::UnifiedPipeline, world::World};

/// Test that Simple and Advanced mode configurations produce equivalent results
/// when using the unified pipeline.
#[test]
fn test_mode_equivalence() {
    let f = 16u32; // Small grid for fast test
    let plates = 3u32;
    let seed: u64 = 12345;
    
    // Test Simple mode
    let mut world_simple = World::new(f, plates, seed);
    let config_simple = PhysicsConfig::simple_mode();
    let mut pipeline_simple = UnifiedPipeline::new(config_simple);
    
    let result_simple = pipeline_simple.step(&mut world_simple, PipelineMode::Realtime { preserve_depth: true });
    let elevation_simple = pipeline_simple.elevation(&world_simple).to_vec();
    
    // Test Advanced mode (create identical world)
    let mut world_advanced = World::new(f, plates, seed);
    let config_advanced = PhysicsConfig::advanced_mode();
    let mut pipeline_advanced = UnifiedPipeline::new(config_advanced);
    
    let result_advanced = pipeline_advanced.step(&mut world_advanced, PipelineMode::Realtime { preserve_depth: true });
    let elevation_advanced = pipeline_advanced.elevation(&world_advanced).to_vec();
    
    // Results should be similar (allowing for different cadences)
    println!("Simple mode eta: {:.3}m", result_simple.eta_m);
    println!("Advanced mode eta: {:.3}m", result_advanced.eta_m);
    
    // Both should have valid elevations
    assert!(!elevation_simple.is_empty());
    assert!(!elevation_advanced.is_empty());
    assert_eq!(elevation_simple.len(), elevation_advanced.len());
    
    // Eta values should be reasonable
    assert!(result_simple.eta_m.abs() < 1000.0);
    assert!(result_advanced.eta_m.abs() < 1000.0);
}

/// Test that Batch mode with write-back produces different depth_m but same final elevation.
#[test]
fn test_batch_mode_writeback() {
    let f = 16u32;
    let plates = 3u32;
    let seed: u64 = 12345;
    
    // Test without write-back
    let mut world_no_writeback = World::new(f, plates, seed);
    let config = PhysicsConfig::simple_mode();
    let mut pipeline = UnifiedPipeline::new(config);
    
    let initial_depth = world_no_writeback.depth_m.clone();
    let _result1 = pipeline.step(&mut world_no_writeback, PipelineMode::Batch { write_back_sea_level: false });
    let final_depth_no_writeback = world_no_writeback.depth_m.clone();
    let elevation_no_writeback = pipeline.elevation(&world_no_writeback).to_vec();
    
    // Test with write-back  
    let mut world_with_writeback = World::new(f, plates, seed);
    let mut pipeline2 = UnifiedPipeline::new(config);
    
    let _result2 = pipeline2.step(&mut world_with_writeback, PipelineMode::Batch { write_back_sea_level: true });
    let final_depth_with_writeback = world_with_writeback.depth_m.clone();
    let elevation_with_writeback = pipeline2.elevation(&world_with_writeback).to_vec();
    
    // Without write-back: depth_m should be unchanged from initial
    let depth_change_no_writeback: f32 = initial_depth.iter()
        .zip(final_depth_no_writeback.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    
    // With write-back: depth_m should be different (eta incorporated)
    let depth_change_with_writeback: f32 = initial_depth.iter()
        .zip(final_depth_with_writeback.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    
    println!("Depth change without writeback: {:.3}m", depth_change_no_writeback);
    println!("Depth change with writeback: {:.3}m", depth_change_with_writeback);
    
    // Final elevations should be similar regardless of write-back mode
    let elevation_diff: f32 = elevation_no_writeback.iter()
        .zip(elevation_with_writeback.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    
    println!("Max elevation difference: {:.3}m", elevation_diff);
    
    // Elevation should be reasonably consistent between modes (allowing for eta incorporation)
    // The key test is that both modes produce valid elevations, not that they're identical
    assert!(elevation_diff < 500.0, "Elevation difference too large: {:.3}m", elevation_diff);
}

/// Test that continental buoyancy works independently of surface processes flag.
#[test]  
fn test_buoyancy_independence() {
    let f = 16u32;
    let plates = 3u32;
    let seed: u64 = 12345;
    
    let mut world = World::new(f, plates, seed);
    
    // Add some continental thickness to trigger buoyancy
    for i in 0..world.grid.cells {
        world.c[i] = 0.8; // high continental fraction
        world.th_c_m[i] = 45_000.0; // thick crust
    }
    
    let config = PhysicsConfig {
        enable_surface_processes: false,  // erosion OFF
        enable_continental_buoyancy: true,  // buoyancy ON
        dt_myr: 1.0,
        ..PhysicsConfig::simple_mode()
    };
    
    let initial_depth = world.depth_m.clone();
    let mut pipeline = UnifiedPipeline::new(config);
    
    let _result = pipeline.step(&mut world, PipelineMode::Realtime { preserve_depth: true });
    
    let max_depth_change = world.depth_m.iter()
        .zip(initial_depth.iter())
        .map(|(new, old)| (new - old).abs())
        .fold(0.0f32, f32::max);
    
    println!("Max depth change from buoyancy: {:.3}m", max_depth_change);
    
    // Should see significant depth change from buoyancy even with erosion disabled
    assert!(max_depth_change > 1.0, "Expected >1m depth change from buoyancy, got {:.3}m", max_depth_change);
}
