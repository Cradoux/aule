# Aulé Architecture Smells

## Duplicate Code Paths

### 1. Dual Pipeline Execution (Critical)
**Files**: `aule/engine/src/world.rs:674` (`step_once`) vs `aule/engine/src/pipeline.rs:225` (`step_full`)

**Symptom**: Two complete physics pipelines with different semantics
- `step_once`: Writes sea-level offset into `depth_m`, renders as `elev = -depth_m`  
- `step_full`: Solves `eta` separately, renders as `elev = eta - depth_m`

**Cause**: Historical evolution - `step_full` was added later for viewer needs without refactoring `step_once`

**Risk**: Renderer inconsistency, coastline flicker, maintenance burden

**Fix**: Unify around single pipeline with configurable output mode

### 2. Continental Advection Algorithms (High)
**Files**: 
- `aule/engine/src/continent.rs:487` (`rigid_advect_c_thc_from`) - push-based with 1-ring mapping
- `aule/engine/src/continent.rs:360` (`advect_c_thc`) - pull-based semi-Lagrangian

**Symptom**: Two different algorithms for the same physical process

**Cause**: Simple mode uses one-shot remap, Advanced mode uses per-step advection

**Risk**: Numerical differences between modes, code duplication

**Fix**: Standardize on one algorithm with appropriate sub-stepping

### 3. CPU vs GPU Flexure Paths (Medium)
**Files**: `aule/engine/src/flexure.rs` vs `aule/engine/src/flexure_gpu.rs`

**Symptom**: GPU flexure only available in Simple mode via `pipeline::step_full`

**Cause**: Advanced mode (`world::step_once`) hardcoded to CPU flexure

**Risk**: Performance disparity between modes

**Fix**: Make flexure backend configurable in both pipelines

## Dead/Unreachable Branches

### 4. Unused Staging Buffers (Low)
**Files**: `aule/viewer/src/main.rs:27` (commented out `ELEV_STAGE`)

**Symptom**: Removed staging buffer but comment remains

**Cause**: Code evolution artifact

**Risk**: Code clarity, potential confusion

**Fix**: Remove dead comments

### 5. Legacy Stepper Module (Low)
**Files**: `aule/engine/src/stepper.rs:1-6`

**Symptom**: Entire module marked as deprecated/retained for tests only

**Cause**: Replaced by unified pipeline but kept for backward compatibility

**Risk**: Maintenance burden, confusion

**Fix**: Remove after migrating remaining test dependencies

## Tight Coupling Issues

### 6. Erosion Flag Gates Unrelated Processes (Critical)
**Files**: `aule/engine/src/pipeline.rs:1721`

**Symptom**: `enable_erosion` flag affects topology freeze detection for unrelated processes

**Cause**: Single flag used to gate multiple unrelated behaviors

**Risk**: Cannot enable buoyancy/isostasy without enabling erosion

**Fix**: Separate flags for each process, remove artificial dependencies

### 7. UI State Directly Controls Physics (Medium)
**Files**: `aule/viewer/src/main.rs:232` (`sp.do_surface` → `enable_erosion`)

**Symptom**: UI overlay state directly maps to physics configuration

**Cause**: No abstraction layer between UI and simulation

**Risk**: UI changes can break physics, tight coupling

**Fix**: Introduce configuration abstraction layer

## Multiple Sources of Truth

### 8. Elevation Buffer Proliferation (High)
**Files**: 
- `aule/engine/src/world.rs:44` (`depth_m`)
- `aule/engine/src/world.rs:75` (`sea.eta_m`) 
- `aule/viewer/src/main.rs:26` (`ELEV_CURR`)
- GPU vertex buffers in globe/raster pipelines

**Symptom**: 4+ different elevation representations with different coordinate systems

**Cause**: Evolutionary complexity, no single source of truth

**Risk**: Inconsistent rendering, synchronization bugs

**Fix**: Single authoritative elevation buffer with consistent coordinate system

### 9. Continental Data Duplication (Medium)
**Files**: 
- `aule/engine/src/world.rs:48-54` (main fields)
- `aule/engine/src/world.rs:49-53` (staging fields)
- `aule/viewer/src/main.rs:40-41` (snapshot fields)

**Symptom**: Continental data replicated across multiple locations

**Cause**: Double buffering and snapshot requirements

**Risk**: Synchronization bugs, memory overhead

**Fix**: Centralized continental data management

## Cross-Layer Knowledge Leaks

### 10. Renderer Knows Pipeline Internals (Medium)
**Files**: `aule/viewer/src/main.rs:2154` (elevation formula: `eta - depth`)

**Symptom**: Renderer hardcodes knowledge of pipeline elevation convention

**Cause**: No abstraction for elevation data access

**Risk**: Renderer breaks when pipeline changes

**Fix**: Elevation accessor API that hides internal representation

### 11. Physics Constants Scattered (Low)
**Files**: 
- `aule/engine/src/lib.rs:82` (`PhysConsts`)
- `aule/engine/src/continent.rs:473-476` (hardcoded densities)
- `aule/engine/src/surface.rs:99-102` (hardcoded caps)

**Symptom**: Physical constants defined in multiple places

**Cause**: No centralized constants management

**Risk**: Inconsistent values, difficult to modify

**Fix**: Centralize all constants in `PhysConsts`

## Hidden Global State

### 12. Unsafe Global Elevation Buffer (High)
**Files**: `aule/viewer/src/main.rs:26` (`static mut ELEV_CURR`)

**Symptom**: Mutable global state accessed via `unsafe`

**Cause**: Cross-thread elevation sharing without proper synchronization

**Risk**: Data races, undefined behavior

**Fix**: Replace with proper thread-safe mechanism (channels, Arc<Mutex<>>)

### 13. Magic Singleton GPU State (Medium)
**Files**: Implicit in GPU buffer management across globe/raster pipelines

**Symptom**: GPU state shared implicitly between rendering paths

**Cause**: No explicit GPU resource management

**Risk**: Resource conflicts, unclear ownership

**Fix**: Explicit GPU resource manager

## Circular Dependencies

### 14. World ↔ Pipeline Circular Reference (Medium)
**Files**: `aule/engine/src/world.rs` imports `pipeline`, `pipeline.rs` imports `world`

**Symptom**: Circular module dependency

**Cause**: `step_once` and `step_full` both operating on same `World` type

**Risk**: Compilation issues, unclear module boundaries

**Fix**: Extract shared types to separate module

## Inconsistent Naming

### 15. Elevation Terminology Chaos (Medium)
**Symptom**: Multiple terms for same concept
- `depth` (positive down)
- `elev` (positive up) 
- `height` (ambiguous)
- `bathymetry` (ocean-specific)
- `hyps` (hypsometry)
- `eta` (sea-level offset)

**Cause**: Evolution without naming standardization

**Risk**: Developer confusion, bugs from sign errors

**Fix**: Standardize on `depth_m` (positive down) and `elevation_m` (positive up)

### 16. Process Flag Naming Inconsistency (Low)
**Files**: `enable_erosion` vs `do_flexure` vs `disable_subduction`

**Symptom**: Inconsistent enable/do/disable prefixes for similar concepts

**Cause**: Different developers, no naming convention

**Risk**: UI confusion, code clarity

**Fix**: Standardize on `enable_*` pattern

## Performance Footguns

### 17. Per-Frame Elevation Recomputation (High)
**Files**: `aule/viewer/src/main.rs:2154` (`curr.iter().map(|&z| world.sea.eta_m - z)`)

**Symptom**: Elevation computed from depth every frame even when unchanged

**Cause**: No caching of elevation computation

**Risk**: Unnecessary CPU overhead

**Fix**: Cache elevation and invalidate on world changes only

### 18. Redundant Continental Advection (Medium)
**Files**: `aule/engine/src/continent.rs:390-448` (distance calculations in inner loop)

**Symptom**: O(N²) distance calculations for continental advection

**Cause**: Brute-force neighbor search without spatial indexing

**Risk**: Poor scaling with grid resolution

**Fix**: Use grid's existing neighbor lists more efficiently

### 19. GPU→CPU→GPU Buffer Copies (Medium)
**Files**: Implicit in elevation snapshot mechanism

**Symptom**: Elevation data copied from GPU to CPU (ELEV_CURR) then back to GPU for rendering

**Cause**: No GPU-to-GPU direct path for elevation updates

**Risk**: Memory bandwidth waste, latency

**Fix**: Keep elevation data on GPU when possible

## Summary by Priority

**Critical (Fix First)**:
- Dual Pipeline Execution (#1)
- Erosion Flag Gates Unrelated Processes (#6)

**High Priority**:
- Continental Advection Algorithms (#2)  
- Elevation Buffer Proliferation (#8)
- Unsafe Global Elevation Buffer (#12)
- Per-Frame Elevation Recomputation (#17)

**Medium Priority**:
- CPU vs GPU Flexure Paths (#3)
- UI State Directly Controls Physics (#7)
- Continental Data Duplication (#9)
- Renderer Knows Pipeline Internals (#10)
- World ↔ Pipeline Circular Reference (#14)
- Inconsistent Elevation Terminology (#15)
- Redundant Continental Advection (#18)
- GPU→CPU→GPU Buffer Copies (#19)

**Low Priority**:
- Unused Staging Buffers (#4)
- Legacy Stepper Module (#5)
- Physics Constants Scattered (#11)
- Magic Singleton GPU State (#13)
- Process Flag Naming Inconsistency (#16)
