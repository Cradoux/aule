use egui::{Slider, Ui};
use egui_plot::{Line, Plot, PlotPoints, Points};

pub struct FlexureUI {
    pub show: bool,
    pub e_pa: f64,
    pub nu: f64,
    pub te_m: f64,
    pub k_npm3: f64,
    pub p_npm: f64,
    pub dx_m: f64,
    pub n: usize,
    // cached
    alpha_km: f64,
    iters: usize,
    resid: f64,
    rms: f64,
    analytic: Vec<[f64; 2]>,
    numeric: Vec<[f64; 2]>,
}

impl Default for FlexureUI {
    fn default() -> Self {
        Self {
            show: false,
            e_pa: 70e9,
            nu: 0.25,
            te_m: 25_000.0,
            k_npm3: 3.0e8,
            p_npm: 1.0e12,
            dx_m: 5_000.0,
            n: 1025,
            alpha_km: 0.0,
            iters: 0,
            resid: 0.0,
            rms: 0.0,
            analytic: Vec::new(),
            numeric: Vec::new(),
        }
    }
}

impl FlexureUI {
    pub fn recompute(&mut self) {
        let d = engine::flexure::d_from_te(self.e_pa, self.nu, self.te_m);
        let a = engine::flexure::alpha(d, self.k_npm3);
        self.alpha_km = a / 1000.0;
        // build RHS and solve
        let q = engine::flexure::line_load_rhs(self.n, self.dx_m, self.p_npm);
        let sol = engine::flexure::solve_flexure_1d(&q, self.dx_m, d, self.k_npm3, 1e-8, 20_000);
        self.iters = sol.iters;
        self.resid = sol.resid;
        // Build plot data
        let i0 = self.n / 2;
        // Analytic dense line over domain
        let mut a_pts: Vec<[f64; 2]> = Vec::with_capacity(self.n * 2);
        for i in 0..self.n {
            let x = (i as f64 - i0 as f64) * self.dx_m;
            let y = engine::flexure::w_line_analytic(x, self.p_npm, d, self.k_npm3);
            a_pts.push([x / 1000.0, y]);
        }
        self.analytic = a_pts;
        // Numeric sparse points
        let stride = (self.n / 256).max(1);
        let mut n_pts: Vec<[f64; 2]> = Vec::new();
        for i in (0..self.n).step_by(stride) {
            let x = (i as f64 - i0 as f64) * self.dx_m;
            n_pts.push([x / 1000.0, sol.w[i]]);
        }
        self.numeric = n_pts;
        // RMS over |x|<=4a
        let xwin = 4.0 * a;
        let imax = ((xwin / self.dx_m).floor() as isize).min(i0 as isize);
        let mut sum2 = 0.0;
        let mut sum2_ref = 0.0;
        let mut count = 0.0;
        for di in -(imax as isize)..=(imax as isize) {
            let i = (i0 as isize + di) as usize;
            let x = (di as f64) * self.dx_m;
            let wa = engine::flexure::w_line_analytic(x, self.p_npm, d, self.k_npm3);
            let wn = sol.w[i];
            sum2 += (wn - wa) * (wn - wa);
            sum2_ref += wa * wa;
            count += 1.0;
        }
        let rms = (sum2 / count).sqrt();
        let ref_rms = (sum2_ref / count).sqrt().max(1e-12);
        self.rms = rms / ref_rms;
    }

    #[allow(dead_code)]
    pub fn ui(&mut self, ui: &mut Ui) {
        if !self.show {
            return;
        }
        ui.label("Flexure (1D)");
        let mut changed = false;
        changed |= ui.add(Slider::new(&mut self.e_pa, 10e9..=200e9).text("E (Pa)")).changed();
        changed |= ui.add(Slider::new(&mut self.nu, 0.15..=0.35).text("nu")).changed();
        changed |= ui.add(Slider::new(&mut self.te_m, 5_000.0..=60_000.0).text("Te (m)")).changed();
        changed |= ui
            .add(Slider::new(&mut self.k_npm3, 1.0e8..=1.0e9).logarithmic(true).text("k (N/m^3)"))
            .changed();
        changed |= ui
            .add(Slider::new(&mut self.p_npm, 1.0e10..=5.0e12).logarithmic(true).text("P (N/m)"))
            .changed();
        changed |= ui.add(Slider::new(&mut self.dx_m, 1000.0..=20_000.0).text("dx (m)")).changed();
        changed |= ui.add(Slider::new(&mut self.n, 257..=4097).text("N")).changed();

        if changed {
            self.recompute();
        }

        ui.label(format!("alpha = {:.1} km", self.alpha_km));
        ui.label(format!("CG iters = {}  resid = {:.2e}", self.iters, self.resid));
        ui.label(format!("RMS (|x|<=4α) = {:.2}%", self.rms * 100.0));

        Plot::new("flexure_plot")
            .view_aspect(2.0)
            .x_axis_label("x (km)")
            .y_axis_label("w (m)")
            .show(ui, |plot_ui| {
                plot_ui.line(Line::new(PlotPoints::from(self.analytic.clone())).name("analytic"));
                plot_ui.points(Points::new(PlotPoints::from(self.numeric.clone())).name("numeric"));
            });
    }
}
