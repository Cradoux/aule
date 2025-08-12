use egui::Color32;
use egui_plot::{Line, PlotPoints, Points};

pub struct AgeDepthPlotParams {
    pub sample_cap: usize,
    pub bin_width_myr: f32,
}

pub struct AgeDepthStats {
    pub rms_m: f32,
    pub n_samples: usize,
    pub age_minmax: (f32, f32),
    pub depth_minmax: (f32, f32),
}

pub struct AgeDepthPlotData {
    pub scatter: Points,
    pub binned: Line,
    pub reference: Line,
    pub stats: AgeDepthStats,
}

fn subsample_stride(n: usize, cap: usize) -> usize {
    (n / cap.max(1)).max(1)
}

pub fn build_age_depth_plot(
    age_myr: &[f32],
    depth_m: &[f32],
    params: AgeDepthPlotParams,
    depth_ref_fn: &dyn Fn(f64) -> f64,
) -> AgeDepthPlotData {
    let n = age_myr.len().min(depth_m.len());
    let stride = subsample_stride(n, params.sample_cap);
    let mut pts: Vec<[f64; 2]> = Vec::with_capacity((n + stride - 1) / stride);
    let (mut amin, mut amax) = (f32::INFINITY, f32::NEG_INFINITY);
    let (mut dmin, mut dmax) = (f32::INFINITY, f32::NEG_INFINITY);
    let mut sq_err_sum: f64 = 0.0;
    let mut count: usize = 0;
    for i in (0..n).step_by(stride) {
        let a = age_myr[i];
        if !a.is_finite() || a <= 0.0 {
            continue;
        }
        let d = depth_m[i];
        if !d.is_finite() {
            continue;
        }
        amin = amin.min(a);
        amax = amax.max(a);
        dmin = dmin.min(d);
        dmax = dmax.max(d);
        let dref = depth_ref_fn(a as f64) as f32;
        let err = (d - dref) as f64;
        sq_err_sum += err * err;
        count += 1;
        pts.push([a as f64, d as f64]);
    }

    let rms = if count > 0 { (sq_err_sum / (count as f64)).sqrt() as f32 } else { 0.0 };
    let scatter = Points::new(PlotPoints::from_iter(pts.iter().copied()))
        .color(Color32::from_rgb(80, 160, 240))
        .radius(1.5)
        .name("samples");

    // Binned means
    let bw = params.bin_width_myr.max(1.0);
    let mut bins: Vec<(f32, f32, u32)> = Vec::new();
    if amin.is_finite() && amax.is_finite() && amin < amax {
        let nb = ((amax - amin) / bw).ceil() as usize;
        bins.resize(nb, (0.0, 0.0, 0));
        for i in (0..n).step_by(stride) {
            let a = age_myr[i];
            if a <= 0.0 || !a.is_finite() {
                continue;
            }
            let d = depth_m[i];
            if !d.is_finite() {
                continue;
            }
            let bi = ((a - amin) / bw).floor() as isize;
            if bi >= 0 && (bi as usize) < nb {
                let (s, c, n) = bins[bi as usize];
                bins[bi as usize] = (s + d, c + a, n + 1);
            }
        }
    }
    let mut bpts: Vec<[f64; 2]> = Vec::new();
    for (sum_d, sum_a, n) in bins.iter().copied() {
        if n == 0 {
            continue;
        }
        let mean_d = sum_d / (n as f32);
        let mean_a = sum_a / (n as f32);
        bpts.push([mean_a as f64, mean_d as f64]);
    }
    let binned = Line::new(PlotPoints::from_iter(bpts.iter().copied()))
        .color(Color32::from_rgb(250, 200, 80))
        .name("binned mean");

    // Reference curve densely sampled
    let mut rpts: Vec<[f64; 2]> = Vec::new();
    if amin.is_finite() && amax.is_finite() && amin < amax {
        let steps = 256;
        for k in 0..=steps {
            let t = k as f64 / steps as f64;
            let a = (1.0 - t) * (amin as f64) + t * (amax as f64);
            let d = depth_ref_fn(a);
            rpts.push([a, d]);
        }
    }
    let reference = Line::new(PlotPoints::from_iter(rpts.iter().copied()))
        .color(Color32::WHITE)
        .name("reference");

    AgeDepthPlotData {
        scatter,
        binned,
        reference,
        stats: AgeDepthStats {
            rms_m: rms,
            n_samples: count,
            age_minmax: (amin, amax),
            depth_minmax: (dmin, dmax),
        },
    }
}
