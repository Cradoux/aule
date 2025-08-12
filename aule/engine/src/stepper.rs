//! World stepper (T-400): update world by one kinematic step.

use crate::{age::depth_from_age, boundaries::Boundaries, ridge, subduction, world::World};

/// Step parameters
pub struct StepParams {
    /// Time step in Myr.
    pub dt_myr: f32,
    /// Opening/closing threshold in m/yr used by boundary classification.
    pub tau_open_m_per_yr: f64,
}

/// Step result stats
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepStats {
    /// Simulation time after the step (Myr).
    pub t_myr: f64,
    /// Step index after increment.
    pub step_idx: u64,
}

/// Execute one step. Order:
/// 1) velocities from plates (kept fixed already in `world.plates`)
/// 2) boundaries classify
/// 3) ridge births
/// 4) age += dt except at ridges
/// 5) subduction edits depth (uses age & boundaries)
/// 6) depth from age curve
pub fn step(world: &mut World, p: &StepParams) -> StepStats {
    // 1) velocities: plates already hold fixed Euler poles; v_en set from plates
    world.v_en.clone_from(&world.plates.vel_en);

    // 2) boundaries classification
    world.boundaries =
        Boundaries::classify(&world.grid, &world.plates.plate_id, &world.v_en, p.tau_open_m_per_yr);

    // 3) ridge births (and fringe disabled for now)
    let mut ages = world.age_myr.clone();
    let _ridge_stats = ridge::apply_ridge(
        &world.grid,
        &world.boundaries,
        &mut ages,
        ridge::RidgeParams { fringe_age_myr: 0.0 },
    );

    // 4) age growth except ridge cells set to 0 this step
    let dt = p.dt_myr;
    for (aw, ar) in world.age_myr.iter_mut().zip(ages.iter()) {
        if *ar == 0.0 {
            *aw = 0.0;
        } else {
            *aw += dt;
        }
    }

    // 5) subduction bathy edits (apply deltas on top of curve; we will overwrite depth below)
    let mut depth_temp = vec![0.0f32; world.grid.cells];
    let _sub = subduction::apply_subduction(
        &world.grid,
        &world.boundaries,
        &world.plates.plate_id,
        &world.age_myr,
        &world.v_en,
        &mut depth_temp,
        subduction::SubductionParams {
            tau_conv_m_per_yr: p.tau_open_m_per_yr,
            trench_half_width_km: 50.0,
            arc_offset_km: 150.0,
            arc_half_width_km: 30.0,
            backarc_width_km: 150.0,
            trench_deepen_m: 3000.0,
            arc_uplift_m: -500.0,
            backarc_uplift_m: -200.0,
            rollback_offset_m: 0.0,
            rollback_rate_km_per_myr: 0.0,
            backarc_extension_mode: false,
            backarc_extension_deepen_m: 600.0,
        },
    );

    // 6) depth from age mapping (overwrites)
    for i in 0..world.grid.cells {
        let mut d = depth_from_age(world.age_myr[i] as f64, 2600.0, 350.0, 0.0) as f32;
        if !d.is_finite() {
            d = 6000.0;
        }
        world.depth_m[i] = d.clamp(0.0, 6000.0);
    }

    // clock
    world.clock.t_myr += dt as f64;
    world.clock.step_idx += 1;
    StepStats { t_myr: world.clock.t_myr, step_idx: world.clock.step_idx }
}
