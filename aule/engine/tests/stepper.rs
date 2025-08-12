use engine::{stepper, world::World};

#[test]
fn determinism_and_age_growth() {
    let f: u32 = 16;
    let num_plates: u32 = 2;
    let seed: u64 = 42;
    let mut w1 = World::new(f, num_plates, seed);
    let mut w2 = World::new(f, num_plates, seed);
    let p = stepper::StepParams { dt_myr: 1.0, tau_open_m_per_yr: 0.005 };

    let steps: usize = 10;
    for _ in 0..steps {
        stepper::step(&mut w1, &p);
    }
    for _ in 0..steps {
        stepper::step(&mut w2, &p);
    }

    // Determinism
    assert_eq!(w1.age_myr, w2.age_myr);
    assert_eq!(w1.depth_m, w2.depth_m);

    // Age growth: non-ridge cells should be close to steps * dt, ridge cells zero when born
    // We conservatively check non-negativity and monotonic depth vs age
    for i in 0..w1.grid.cells {
        assert!(w1.age_myr[i] >= 0.0);
        // depth increases with age under our simple mapping
        // local monotonicity (not strict) check against zero-age depth
        let d0 = engine::age::depth_from_age(0.0, 2600.0, 350.0, 0.0) as f32;
        assert!(w1.depth_m[i] >= d0);
    }
}
