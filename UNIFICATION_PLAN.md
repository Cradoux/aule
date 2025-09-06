# Aulé Unification Plan

## Overview

This plan unifies Simple and Advanced modes around a single pipeline architecture while removing architectural smells and improving maintainability. The approach is incremental with each step validated by tests.

## Phase 1: Foundation (Week 1)

### PR1.1: Standardize Elevation Terminology
**Scope**: Rename variables and functions for consistent naming
**Files**: `world.rs`, `pipeline.rs`, viewer files
**Changes**:
- Standardize on `depth_m` (positive down, primary storage)
- Standardize on `elevation_m` (positive up, computed as `-depth_m + eta_m`)
- Remove ambiguous terms: `elev`, `height`, `hyps`, `bathym`
**Risk**: Low - pure refactoring
**Test Gate**: All existing tests pass with renamed variables

### PR1.2: Extract Shared Configuration Types
**Scope**: Break circular dependency between `world.rs` and `pipeline.rs`
**Files**: New `config.rs`, `world.rs`, `pipeline.rs`
**Changes**:
- Create `config.rs` with `PhysicsConfig`, `RenderConfig` types
- Move `PipelineCfg` and `StepParams` to config module
- Remove circular imports
**Risk**: Low - pure refactoring
**Test Gate**: Compilation succeeds, no behavior change

### PR1.3: Centralize Physical Constants
**Scope**: Single source of truth for all physical constants
**Files**: `lib.rs`, `continent.rs`, `surface.rs`, others
**Changes**:
- Extend `PhysConsts` with all hardcoded values
- Replace scattered constants with `PhysConsts` references
- Add validation for constant ranges
**Risk**: Low - constants should not change behavior
**Test Gate**: Numerical outputs identical before/after

## Phase 2: Pipeline Unification (Week 2)

### PR2.1: Create Single Pipeline Interface
**Scope**: Abstract pipeline execution behind unified API
**Files**: New `unified_pipeline.rs`, `world.rs`, `pipeline.rs`
**Changes**:
```rust
pub struct UnifiedPipeline {
    config: PhysicsConfig,
}

pub enum PipelineMode {
    Realtime { preserve_depth: bool },
    Batch { write_back_sea_level: bool },
}

impl UnifiedPipeline {
    pub fn step(&mut self, world: &mut World, mode: PipelineMode) -> StepResult;
}
```
**Risk**: Medium - new abstraction layer
**Test Gate**: Both modes produce identical results to current implementation

### PR2.2: Unify Sea-Level Handling
**Scope**: Single sea-level calculation with configurable output
**Files**: `unified_pipeline.rs`, `isostasy.rs`
**Changes**:
- Always solve `eta_m` using isostasy logic
- Mode flag controls whether to write back to `depth_m` or keep separate
- Renderer always uses `elevation = -depth_m + eta_m`
**Risk**: High - changes core elevation semantics
**Test Gate**: Visual comparison shows identical coastlines in both modes

### PR2.3: Standardize Continental Advection
**Scope**: Single advection algorithm with sub-stepping
**Files**: `continent.rs`
**Changes**:
- Keep `rigid_advect_c_thc_from` (more robust) as canonical implementation
- Remove `advect_c_thc` (less stable)
- Add sub-stepping parameter for stability
- Use same algorithm in both modes
**Risk**: Medium - changes numerical behavior
**Test Gate**: Continental drift rates within 5% of current Advanced mode

## Phase 3: Buffer Management (Week 2-3)

### PR3.1: Eliminate Unsafe Global State
**Scope**: Replace `static mut ELEV_CURR` with thread-safe mechanism
**Files**: `viewer/src/main.rs`
**Changes**:
- Replace with `Arc<RwLock<Vec<f32>>>` or channel-based system
- Update all access points to use safe API
- Add proper error handling for lock contention
**Risk**: Medium - threading behavior change
**Test Gate**: No crashes under concurrent access, same visual output

### PR3.2: Centralize Elevation Buffer Management
**Scope**: Single source of truth for elevation data
**Files**: New `elevation_manager.rs`, viewer files
**Changes**:
```rust
pub struct ElevationManager {
    depth_m: Vec<f32>,
    eta_m: f32,
    cached_elevation: Option<Vec<f32>>,
    cache_valid: bool,
}

impl ElevationManager {
    pub fn get_elevation(&mut self) -> &[f32];
    pub fn invalidate(&mut self);
    pub fn update_depth(&mut self, new_depth: Vec<f32>);
    pub fn update_eta(&mut self, new_eta: f32);
}
```
**Risk**: Medium - centralized state management
**Test Gate**: Elevation values identical, no performance regression

### PR3.3: GPU Buffer Consistency  
**Scope**: Eliminate redundant CPU↔GPU transfers
**Files**: `raster_gpu.rs`, `globe/pipeline.rs`
**Changes**:
- GPU buffers updated only when `ElevationManager` invalidates cache
- Direct GPU-to-GPU copies where possible
- Reduce CPU readback operations
**Risk**: Low - performance optimization
**Test Gate**: Visual output identical, frame time improved

## Phase 4: Process Decoupling (Week 3)

### PR4.1: Separate Process Enable Flags
**Scope**: Independent control of each physics process
**Files**: `config.rs`, `unified_pipeline.rs`, viewer UI
**Changes**:
```rust
pub struct PhysicsConfig {
    pub enable_rigid_motion: bool,
    pub enable_subduction: bool,
    pub enable_transforms: bool,
    pub enable_flexure: bool,
    pub enable_surface_processes: bool,
    pub enable_isostasy: bool,
    pub enable_continental_buoyancy: bool,  // NEW: separate from erosion
}
```
**Risk**: Low - better separation of concerns
**Test Gate**: Can enable buoyancy without enabling erosion

### PR4.2: Flexible Flexure Backend
**Scope**: GPU/CPU flexure selection in both modes
**Files**: `flexure_manager.rs`, `unified_pipeline.rs`
**Changes**:
```rust
pub enum FlexureBackend {
    CpuWinkler,
    GpuMultigrid { levels: u32, cycles: u32 },
}

pub struct FlexureManager {
    backend: FlexureBackend,
    gpu_available: bool,
}
```
**Risk**: Medium - backend abstraction
**Test Gate**: Both backends produce equivalent results within 5%

### PR4.3: Configurable Process Cadences
**Scope**: Unified cadence system for all processes
**Files**: `config.rs`, `unified_pipeline.rs`
**Changes**:
- Single `CadenceConfig` struct with per-process intervals
- Replace scattered cadence logic with centralized scheduler
- Support sub-stepping for stability
**Risk**: Low - cleaner configuration
**Test Gate**: Same effective update rates as before

## Phase 5: Mode Convergence (Week 3)

### PR5.1: Mode as Configuration
**Scope**: Simple/Advanced as parameter sets, not code paths
**Files**: `config.rs`, viewer UI
**Changes**:
```rust
impl PhysicsConfig {
    pub fn simple_mode() -> Self {
        Self {
            cadence_all: 1,  // every step
            flexure_backend: FlexureBackend::GpuMultigrid { levels: 3, cycles: 2 },
            preserve_depth_buffer: true,  // don't write eta back to depth
            // ...
        }
    }
    
    pub fn advanced_mode() -> Self {
        Self {
            cadence_transforms: 2,
            cadence_subduction: 4,
            flexure_backend: FlexureBackend::CpuWinkler,
            preserve_depth_buffer: false,  // write eta back to depth
            // ...
        }
    }
}
```
**Risk**: Low - configuration change only
**Test Gate**: Both modes produce visually identical results

### PR5.2: Remove Legacy Pipeline
**Scope**: Delete `world::step_once`, use only unified pipeline
**Files**: `world.rs`, tests
**Changes**:
- Mark `step_once` as deprecated
- Update all call sites to use `UnifiedPipeline::step`
- Migrate test dependencies
- Remove after migration complete
**Risk**: High - removes major code path
**Test Gate**: All tests pass with unified pipeline

## Phase 6: Performance & Polish (Week 4)

### PR6.1: Elevation Caching Optimization
**Scope**: Eliminate per-frame elevation recomputation
**Files**: `elevation_manager.rs`, viewer rendering
**Changes**:
- Cache computed elevation until invalidated
- Track dirty flags for depth/eta changes
- Batch updates to minimize cache invalidation
**Risk**: Low - performance optimization
**Test Gate**: 20%+ frame time improvement, visual output identical

### PR6.2: Continental Advection Optimization  
**Scope**: Use grid neighbor lists efficiently
**Files**: `continent.rs`
**Changes**:
- Replace O(N²) distance search with grid neighbor traversal
- Pre-compute rotation matrices per plate
- Vectorize distance calculations where possible
**Risk**: Medium - algorithm optimization
**Test Gate**: 50%+ advection time improvement, results within 1% tolerance

### PR6.3: Memory Layout Optimization
**Scope**: Structure-of-Arrays (SoA) consistency
**Files**: `world.rs`, field access patterns
**Changes**:
- Ensure all per-cell data uses SoA layout
- Align buffer sizes for SIMD operations
- Reduce memory fragmentation
**Risk**: Low - memory layout change
**Test Gate**: Memory usage reduced by 15%, no performance regression

## Core API Specifications

### Unified Pipeline API
```rust
pub struct UnifiedPipeline {
    config: PhysicsConfig,
    elevation_mgr: ElevationManager,
    flexure_mgr: FlexureManager,
}

impl UnifiedPipeline {
    /// Execute one physics step with configurable output mode
    pub fn step(&mut self, world: &mut World, mode: PipelineMode) -> StepResult {
        // 1. Boundary classification
        // 2. Kinematics/advection  
        // 3. Ridge birth
        // 4. Subduction/transforms (if enabled)
        // 5. Continental buoyancy (if enabled)
        // 6. Flexure (if enabled, backend configurable)
        // 7. Surface processes (if enabled)
        // 8. Sea-level regulation (always, output mode configurable)
        // 9. Elevation cache invalidation
    }
}
```

### Elevation Management API
```rust
pub trait ElevationProvider {
    /// Get current elevation field (positive up)
    fn elevation(&mut self) -> &[f32];
    /// Check if elevation data has changed since last access
    fn is_dirty(&self) -> bool;
    /// Mark elevation as needing recomputation
    fn invalidate(&mut self);
}
```

### Continental Remapping API  
```rust
pub fn remap_continental_fields(
    grid: &Grid,
    plates: &Plates,
    dt_years: f64,
    src_c: &[f32],
    src_th_c: &[f32],
    dst_c: &mut [f32],
    dst_th_c: &mut [f32],
) -> RemapStats;
```

### Flexure Backend API
```rust
pub trait FlexureBackend {
    fn solve_deflection(
        &mut self,
        grid: &Grid,
        loads: &[f32],
        te_field: &[f32],
        deflection: &mut [f32],
    ) -> FlexureStats;
}
```

## Migration Timeline

| Week | Phase | Deliverables | Risk | Test Gate |
|------|-------|--------------|------|-----------|
| 1 | Foundation | Terminology, Config, Constants | Low | No behavior change |
| 2 | Pipeline | Unified API, Sea-level, Advection | High | Visual parity |
| 2-3 | Buffers | Safe globals, Centralized elevation, GPU consistency | Medium | Performance improvement |
| 3 | Decoupling | Separate flags, Flexible flexure, Cadences | Medium | Independent process control |
| 3 | Convergence | Mode as config, Remove legacy | High | Single pipeline works |
| 4 | Polish | Caching, Optimization, Memory | Low | Performance gains |

## Success Criteria

### Functional Requirements
- [ ] Single pipeline supports both Simple and Advanced behavior
- [ ] Continental drift aligns with plate motion (<0.25° error)
- [ ] Buoyancy works independently of erosion flag
- [ ] Renderer reads consistent elevation buffer in both modes
- [ ] Sea-level freeze affects only flooding, not elevation values

### Performance Requirements  
- [ ] Frame time maintained or improved
- [ ] Memory usage reduced by 10%+
- [ ] Continental advection 50%+ faster
- [ ] GPU utilization available in both modes

### Code Quality Requirements
- [ ] No `unsafe` code in viewer
- [ ] No circular module dependencies  
- [ ] Single source of truth for elevation
- [ ] Consistent naming throughout
- [ ] All processes independently configurable

## Risk Mitigation

### High-Risk Changes
1. **Sea-level unification**: Extensive visual testing, A/B comparison tool
2. **Pipeline removal**: Gradual deprecation, parallel testing period
3. **Elevation semantics**: Automated visual diff testing

### Rollback Strategy
- Each PR is independently revertible
- Feature flags for new pipeline during transition
- Parallel old/new pipeline support during migration
- Automated regression detection

### Testing Strategy
- Unit tests for each core API
- Integration tests for mode parity
- Visual regression tests for coastlines
- Performance benchmarks for each optimization
- Memory leak detection in long-running tests
