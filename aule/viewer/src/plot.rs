use egui::Color32;
use egui_plot::{Line, PlotPoints, Points};
use std::fs;

#[allow(dead_code)]
pub struct AgeDepthPlotParams {
    pub sample_cap: usize,
    pub bin_width_myr: f32,
}

#[allow(dead_code)]
pub struct AgeDepthStats {
    pub rms_m: f32,
    pub n_samples: usize,
    pub age_minmax: (f32, f32),
    pub depth_minmax: (f32, f32),
}

#[allow(dead_code)]
pub struct AgeDepthPlotData {
    pub scatter: Points,
    pub binned: Line,
    pub reference: Line,
    pub stats: AgeDepthStats,
}

#[allow(dead_code)]
pub struct Dl1Series {
    pub t_myr: Vec<f64>,
    pub land_pct: Vec<f64>,
    pub cap_comp: Vec<f64>,
    pub cap_thc: Vec<f64>,
}

#[allow(dead_code)]
pub fn read_dl1_series_csv(path: &str) -> Option<Dl1Series> {
    let s = fs::read_to_string(path).ok()?;
    let mut t = Vec::new();
    let mut land = Vec::new();
    let mut capc = Vec::new();
    let mut capth = Vec::new();
    for line in s.lines() {
        if line.starts_with("t_myr") {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 9 {
            continue;
        }
        let tm = cols[0].parse::<f64>().ok()?;
        let pr_r = cols[1].parse::<f64>().unwrap_or(0.0);
        let pr_s = cols[2].parse::<f64>().unwrap_or(0.0);
        let pr_t = cols[3].parse::<f64>().unwrap_or(0.0);
        let _oceanized = cols[4].parse::<f64>().unwrap_or(0.0);
        let cc = cols[5].parse::<f64>().unwrap_or(0.0);
        let ct = cols[6].parse::<f64>().unwrap_or(0.0);
        // land% approx from ridge/subd/transform shares is not available; compute from world in UI if needed.
        // For CSV we just plot caps by time; land fraction can be computed separately and appended if desired.
        t.push(tm);
        land.push((pr_r + pr_s + pr_t) * 0.0); // placeholder 0, UI can provide live land%
        capc.push(cc);
        capth.push(ct);
    }
    Some(Dl1Series { t_myr: t, land_pct: land, cap_comp: capc, cap_thc: capth })
}

#[allow(dead_code)]
fn subsample_stride(n: usize, cap: usize) -> usize {
    (n / cap.max(1)).max(1)
}

#[allow(dead_code)]
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
