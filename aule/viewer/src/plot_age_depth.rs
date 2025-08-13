use egui::Ui;
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelKind {
    HalfSpace,
    Plate,
}

pub struct AgeDepthUIState {
    pub show: bool,
    pub model: ModelKind,
    pub d0_m: f64,
    pub k_m2_s: f64,
    pub t_plate_myr: f64,
    pub z_plate_m: f64,
    pub sample_cap: usize,
    pub use_final_depth: bool,
    // cache
    last_hash: u64,
    scatter_pts: Option<Vec<[f64; 2]>>,
    reference_pts: Option<Vec<[f64; 2]>>,
    binned_pts: Option<Vec<[f64; 2]>>,
    rms_m: f32,
    n_samples: usize,
    age_minmax: (f32, f32),
    depth_minmax: (f32, f32),
}

impl Default for AgeDepthUIState {
    fn default() -> Self {
        Self {
            show: false,
            model: ModelKind::HalfSpace,
            d0_m: 2600.0,
            k_m2_s: 1.0e-6,
            t_plate_myr: 120.0,
            z_plate_m: 6000.0,
            sample_cap: 10_000,
            use_final_depth: true,
            last_hash: 0,
            scatter_pts: None,
            reference_pts: None,
            binned_pts: None,
            rms_m: 0.0,
            n_samples: 0,
            age_minmax: (0.0, 0.0),
            depth_minmax: (0.0, 0.0),
        }
    }
}

fn compute_hash(params: (&AgeDepthUIState, usize, usize)) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    params.hash(&mut h);
    h.finish()
}

impl std::hash::Hash for AgeDepthUIState {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.show.hash(state);
        (match self.model {
            ModelKind::HalfSpace => 0u8,
            ModelKind::Plate => 1,
        })
        .hash(state);
        // Quantize floats for stability
        ((self.d0_m * 10.0).round() as i64).hash(state);
        ((self.k_m2_s * 1.0e9).round() as i64).hash(state);
        ((self.t_plate_myr * 10.0).round() as i64).hash(state);
        ((self.z_plate_m * 1.0).round() as i64).hash(state);
        self.sample_cap.hash(state);
        self.use_final_depth.hash(state);
    }
}

pub fn ui(
    st: &mut AgeDepthUIState,
    ui: &mut Ui,
    grid: &engine::grid::Grid,
    age_myr: &[f32],
    depth_pre_m: &[f32],
    depth_final_m: &[f32],
) {
    ui.collapsing("Age–Depth (A)", |ui| {
        let mut changed = false;
        ui.horizontal(|ui| {
            changed |= ui.radio_value(&mut st.model, ModelKind::HalfSpace, "Half-space").changed();
            changed |= ui.radio_value(&mut st.model, ModelKind::Plate, "Plate").changed();
        });
        ui.horizontal(|ui| {
            changed |=
                ui.add(egui::Slider::new(&mut st.d0_m, 2000.0..=3000.0).text("d0 (m)")).changed();
            changed |= ui
                .add(
                    egui::Slider::new(&mut st.k_m2_s, 5.0e-7..=2.0e-6)
                        .logarithmic(true)
                        .text("k (m^2/s)"),
                )
                .changed();
            if st.model == ModelKind::Plate {
                changed |= ui
                    .add(egui::Slider::new(&mut st.t_plate_myr, 60.0..=200.0).text("t_plate (Myr)"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut st.z_plate_m, 4000.0..=7000.0).text("z_plate (m)"))
                    .changed();
            }
        });
        ui.horizontal(|ui| {
            changed |=
                ui.checkbox(&mut st.use_final_depth, "Sampling: Final (L+continents)").changed();
            changed |= ui
                .add(
                    egui::Slider::new(&mut st.sample_cap, 1_000..=50_000)
                        .text("Samples")
                        .step_by(1000.0),
                )
                .changed();
        });

        // Recompute if params changed
        let hash = compute_hash((st, grid.cells, st.sample_cap));
        if changed || hash != st.last_hash || st.scatter_pts.is_none() || st.reference_pts.is_none()
        {
            st.last_hash = hash;
            let depth_src = if st.use_final_depth { depth_final_m } else { depth_pre_m };
            let n = grid.cells.min(depth_src.len()).min(age_myr.len());
            let stride = (n / st.sample_cap.max(1)).max(1);
            // Sample points deterministically
            let mut pts: Vec<[f64; 2]> = Vec::with_capacity((n + stride - 1) / stride);
            let (mut amin, mut amax) = (f32::INFINITY, f32::NEG_INFINITY);
            let (mut dmin, mut dmax) = (f32::INFINITY, f32::NEG_INFINITY);
            let mut sq_err_sum: f64 = 0.0;
            let mut count: usize = 0;
            for i in (0..n).step_by(stride) {
                let a = age_myr[i];
                let d = depth_src[i];
                if !a.is_finite() || !d.is_finite() || a <= 0.0 {
                    continue;
                }
                amin = amin.min(a);
                amax = amax.max(a);
                dmin = dmin.min(d);
                dmax = dmax.max(d);
                let d_ref = match st.model {
                    ModelKind::HalfSpace => {
                        engine::age::depth_from_age_hsc(a as f64, st.d0_m, st.k_m2_s)
                    }
                    ModelKind::Plate => engine::age::depth_from_age_plate(
                        a as f64,
                        st.d0_m,
                        st.t_plate_myr,
                        st.z_plate_m,
                        st.k_m2_s,
                    ),
                } as f32;
                let e = d - d_ref;
                sq_err_sum += (e as f64) * (e as f64);
                count += 1;
                pts.push([a as f64, d as f64]);
            }
            st.rms_m = if count > 0 { (sq_err_sum / (count as f64)).sqrt() as f32 } else { 0.0 };
            st.n_samples = count;
            st.age_minmax = (amin, amax);
            st.depth_minmax = (dmin, dmax);
            st.scatter_pts = Some(pts);

            // Build analytic reference curve on [0, 120] Myr
            let mut rpts: Vec<[f64; 2]> = Vec::new();
            let a_lo = 0.0_f64;
            let a_hi = 120.0_f64;
            for k in 0..=480 {
                let t = k as f64 / 480.0;
                let a = (1.0 - t) * a_lo + t * a_hi;
                let d = match st.model {
                    ModelKind::HalfSpace => engine::age::depth_from_age_hsc(a, st.d0_m, st.k_m2_s),
                    ModelKind::Plate => engine::age::depth_from_age_plate(
                        a,
                        st.d0_m,
                        st.t_plate_myr,
                        st.z_plate_m,
                        st.k_m2_s,
                    ),
                };
                rpts.push([a, d]);
            }
            st.reference_pts = Some(rpts);

            // Simple binned mean over samples in [0, 120]
            let mut bins: Vec<(f64, f64, u32)> = vec![(0.0, 0.0, 0); 120 / 2];
            for i in (0..n).step_by(stride) {
                let a = age_myr[i] as f64;
                let d = depth_src[i] as f64;
                if !(a.is_finite() && d.is_finite()) || a <= 0.0 || a > 120.0 {
                    continue;
                }
                let bi = (a / 2.0).floor() as isize;
                if bi >= 0 && (bi as usize) < bins.len() {
                    let (sd, sa, c) = bins[bi as usize];
                    bins[bi as usize] = (sd + d, sa + a, c + 1);
                }
            }
            let mut bpts: Vec<[f64; 2]> = Vec::new();
            for (sd, sa, c) in bins {
                if c == 0 {
                    continue;
                }
                bpts.push([sa / (c as f64), sd / (c as f64)]);
            }
            st.binned_pts = Some(bpts);
        }

        // Draw
        Plot::new("age_depth_validation")
            .legend(Legend::default())
            .x_axis_label("Age (Myr)")
            .y_axis_label("Depth (m)")
            .show(ui, |plot_ui| {
                if let Some(r) = &st.reference_pts {
                    plot_ui
                        .line(Line::new(PlotPoints::from_iter(r.iter().copied())).name("analytic"));
                }
                if let Some(p) = &st.scatter_pts {
                    plot_ui.points(
                        Points::new(PlotPoints::from_iter(p.iter().copied())).name("samples"),
                    );
                }
                if let Some(b) = &st.binned_pts {
                    plot_ui
                        .line(Line::new(PlotPoints::from_iter(b.iter().copied())).name("binned"));
                }
            });
        ui.label(format!(
            "RMS: {:.1} m  N={}  Age[min/max]={:.2}/{:.2} Myr  Depth[min/max]={:.0}/{:.0} m",
            st.rms_m,
            st.n_samples,
            st.age_minmax.0,
            st.age_minmax.1,
            st.depth_minmax.0,
            st.depth_minmax.1
        ));
    });
}
