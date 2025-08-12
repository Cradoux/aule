//! WGSL multigrid scaffold for thin-plate + Winkler on tiled atlas (cartesian per-tile).
//! This initial version exposes the A-operator and a minimal V-cycle shell.
#![allow(missing_docs)]

use crate::gpu::GpuContext;
use wgpu::util::DeviceExt;

/// Tile dimensions (including halo depth `halo` on each side)
#[derive(Clone, Copy, Debug)]
pub struct TileDims {
    /// Interior width (without halos)
    pub width: u32,
    /// Interior height (without halos)
    pub height: u32,
    /// Halo cells on each side
    pub halo: u32,
    /// Grid spacing (m)
    pub dx: f32,
}

impl TileDims {
    /// Total width including halos
    pub fn width_total(&self) -> u32 {
        self.width + 2 * self.halo
    }
    /// Total height including halos
    pub fn height_total(&self) -> u32 {
        self.height + 2 * self.halo
    }
}

/// Parameters for flexure solve
#[derive(Clone, Copy, Debug)]
pub struct FlexParams {
    /// Plate rigidity D (Pa·m^3)
    pub d: f32,
    /// Winkler stiffness k (N/m^3)
    pub k: f32,
    /// Weighted-Jacobi relaxation factor
    pub wj_omega: f32,
    /// Pre-smoothing sweeps
    pub nu1: u32,
    /// Post-smoothing sweeps
    pub nu2: u32,
    /// Number of V-cycle levels
    pub levels: u32,
}

/// Residual norms
#[derive(Clone, Copy, Debug)]
pub struct FlexStats {
    /// L2 residual before cycle
    pub res_in: f64,
    /// L2 residual after cycle
    pub res_out: f64,
}

/// Simple buffer wrapper
pub struct GpuTex {
    pub buf: wgpu::Buffer,
    pub len: usize,
}

impl GpuTex {
    pub fn new_storage(ctx: &GpuContext, len: usize, label: &str) -> Self {
        let buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: (len * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        Self { buf, len }
    }
}

/// Scratch buffers for one level (kept simple for now)
pub struct FlexScratch {
    pub r: GpuTex,
    pub tmp: GpuTex,
}

impl FlexScratch {
    pub fn new(ctx: &GpuContext, n: usize) -> Self {
        Self {
            r: GpuTex::new_storage(ctx, n, "flex.r"),
            tmp: GpuTex::new_storage(ctx, n, "flex.tmp"),
        }
    }
}

/// GPU pipelines for flexure operator and simple components
pub struct FlexGpu {
    pipeline_apply_a: wgpu::ComputePipeline,
    bind_apply_a: wgpu::BindGroupLayout,
}

impl FlexGpu {
    /// Create pipelines
    pub fn new(ctx: &GpuContext) -> Self {
        let shader_src =
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../shaders/flexure.wgsl"));
        let module = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("flexure.wgsl"),
            source: wgpu::ShaderSource::Wgsl(shader_src.into()),
        });
        let bind_apply_a = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bind.apply_a"),
            entries: &[
                // params uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(32),
                    },
                    count: None,
                },
                // in w
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // out Aw
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });
        let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pipe.flexure"),
            bind_group_layouts: &[&bind_apply_a],
            push_constant_ranges: &[],
        });
        let pipeline_apply_a =
            ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("apply_A"),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: "apply_A",
            });
        Self { pipeline_apply_a, bind_apply_a }
    }

    fn make_params_buf(&self, ctx: &GpuContext, dims: TileDims, p: &FlexParams) -> wgpu::Buffer {
        #[repr(C)]
        #[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
        struct Params {
            w_tot: u32,
            h_tot: u32,
            halo: u32,
            _pad: u32,
            dx: f32,
            d: f32,
            k: f32,
            _pad2: f32,
        }
        let params = Params {
            w_tot: dims.width_total(),
            h_tot: dims.height_total(),
            halo: dims.halo,
            _pad: 0,
            dx: dims.dx,
            d: p.d,
            k: p.k,
            _pad2: 0.0,
        };
        ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("params"),
            contents: bytemuck::cast_slice(std::slice::from_ref(&params)),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        })
    }

    /// Apply A into out buffer: out = D*(D4x+D4y) + k*w
    pub fn dispatch_apply_a(
        &self,
        ctx: &GpuContext,
        dims: TileDims,
        w_in: &GpuTex,
        out: &GpuTex,
        p: &FlexParams,
    ) {
        let params_buf = self.make_params_buf(ctx, dims, p);
        let bind = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bg.apply_a"),
            layout: &self.bind_apply_a,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: w_in.buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: out.buf.as_entire_binding() },
            ],
        });
        let mut encoder = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("enc.apply_a") });
        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("cpass.apply_a"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&self.pipeline_apply_a);
            cpass.set_bind_group(0, &bind, &[]);
            let gx = ((dims.width) + 7) / 8;
            let gy = ((dims.height) + 7) / 8;
            cpass.dispatch_workgroups(gx, gy, 1);
        }
        ctx.queue.submit(Some(encoder.finish()));
    }

    /// Minimal V-cycle scaffold: compute residual before and after one apply (placeholder smoother)
    pub fn v_cycle(
        &mut self,
        ctx: &GpuContext,
        dims: TileDims,
        w: &mut GpuTex,
        f: &GpuTex,
        tmp: &mut FlexScratch,
        p: &FlexParams,
    ) -> FlexStats {
        // Compute r0 = f - A w
        self.dispatch_apply_a(ctx, dims, w, &tmp.tmp, p);
        // r = f - tmp
        // Do reduction on CPU: map both
        let n = (dims.width_total() * dims.height_total()) as usize;
        let slice_tmp = blocking_read(&ctx.device, &ctx.queue, &tmp.tmp.buf, n);
        let slice_f = blocking_read(&ctx.device, &ctx.queue, &f.buf, n);
        let res_in = l2_residual(&slice_f, &slice_tmp);
        // Dummy smoother: one weighted-Jacobi like step in place: w += omega*(f - A w)/k
        let omega = if p.k > 0.0 { p.wj_omega.clamp(0.0, 1.0) } else { 0.7 };
        let mut w_host = blocking_read(&ctx.device, &ctx.queue, &w.buf, n);
        let k = p.k as f64;
        for i in 0..n {
            let r_i = (slice_f[i] as f64) - (slice_tmp[i] as f64);
            if k > 0.0 {
                w_host[i] = (w_host[i] as f64 + (omega as f64) * r_i / k) as f32;
            }
        }
        ctx.queue.write_buffer(&w.buf, 0, bytemuck::cast_slice(&w_host));
        // Recompute residual
        self.dispatch_apply_a(ctx, dims, w, &tmp.tmp, p);
        let slice_tmp2 = blocking_read(&ctx.device, &ctx.queue, &tmp.tmp.buf, n);
        let res_out = l2_residual(&slice_f, &slice_tmp2);
        FlexStats { res_in, res_out }
    }
}

fn l2_residual(f: &[f32], aw: &[f32]) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..f.len() {
        let r = (f[i] as f64) - (aw[i] as f64);
        sum += r * r;
    }
    sum.sqrt()
}

fn blocking_read(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    buf: &wgpu::Buffer,
    n: usize,
) -> Vec<f32> {
    // Copy to a temporary read buffer then map
    let size = (n * std::mem::size_of::<f32>()) as u64;
    let read_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("readback"),
        size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    let mut encoder =
        device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("enc.read") });
    encoder.copy_buffer_to_buffer(buf, 0, &read_buf, 0, size);
    queue.submit(Some(encoder.finish()));
    // Map
    read_buf.slice(..).map_async(wgpu::MapMode::Read, |_| {});
    device.poll(wgpu::Maintain::Wait);
    let data = read_buf.slice(..).get_mapped_range();
    let mut out = vec![0.0f32; n];
    out.copy_from_slice(bytemuck::cast_slice(&data));
    drop(data);
    read_buf.unmap();
    out
}

// ---------- CPU manufactured solution helpers ----------

/// Build w_true and f = D*(D4x+D4y)w + k w on a cartesian grid with halo=2, Dirichlet outside
pub fn build_manufactured_rhs(
    width: usize,
    height: usize,
    halo: usize,
    dx: f32,
    d: f32,
    k: f32,
) -> (Vec<f32>, Vec<f32>) {
    let wt = width + 2 * halo;
    let ht = height + 2 * halo;
    let mut w = vec![0.0f32; wt * ht];
    let mut f = vec![0.0f32; wt * ht];
    let l = (width as f32) * dx;
    for j in 0..height {
        for i in 0..width {
            let x = (i as f32) * dx;
            let y = (j as f32) * dx;
            let val = (2.0 * std::f32::consts::PI * x / l).sin()
                * (2.0 * std::f32::consts::PI * y / l).sin();
            let idx = (j + halo) * wt + (i + halo);
            w[idx] = val;
        }
    }
    // CPU apply A
    let inv_dx4 = 1.0f32 / (dx * dx * dx * dx);
    for j in halo..(height + halo) {
        for i in halo..(width + halo) {
            let idx = j * wt + i;
            let w_im2 = w[j * wt + (i - 2)];
            let w_im1 = w[j * wt + (i - 1)];
            let w_ip1 = w[j * wt + (i + 1)];
            let w_ip2 = w[j * wt + (i + 2)];
            let d4x = (w_im2 - 4.0 * w_im1 + 6.0 * w[idx] - 4.0 * w_ip1 + w_ip2) * inv_dx4;
            let w_jm2 = w[(j - 2) * wt + i];
            let w_jm1 = w[(j - 1) * wt + i];
            let w_jp1 = w[(j + 1) * wt + i];
            let w_jp2 = w[(j + 2) * wt + i];
            let d4y = (w_jm2 - 4.0 * w_jm1 + 6.0 * w[idx] - 4.0 * w_jp1 + w_jp2) * inv_dx4;
            f[idx] = d * (d4x + d4y) + k * w[idx];
        }
    }
    (w, f)
}
