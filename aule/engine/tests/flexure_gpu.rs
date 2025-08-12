use engine::flexure_gpu::{
    build_manufactured_rhs, FlexGpu, FlexParams, FlexScratch, GpuTex, TileDims,
};

#[test]
fn cpu_operator_matches_wgsl_on_small_patch() {
    // 16x16 interior + halo=2
    let width = 16usize;
    let height = 16usize;
    let halo = 2usize;
    let dx = 1.0f32;
    let d = 2.0f32;
    let k = 3.0f32;
    let (w_host, f_host) = build_manufactured_rhs(width, height, halo, dx, d, k);
    let wt = width + 2 * halo;
    let ht = height + 2 * halo;
    let n = wt * ht;

    let ctx = pollster::block_on(engine::gpu::GpuContext::new());
    let flex = FlexGpu::new(&ctx);
    let dims = TileDims { width: width as u32, height: height as u32, halo: halo as u32, dx };
    let p = FlexParams { d, k, wj_omega: 0.7, nu1: 2, nu2: 2, levels: 1 };

    let w = {
        let t = GpuTex::new_storage(&ctx, n, "w");
        ctx.queue.write_buffer(&t.buf, 0, bytemuck::cast_slice(&w_host));
        t
    };
    let out = GpuTex::new_storage(&ctx, n, "out");

    flex.dispatch_apply_a(&ctx, dims, &w, &out, &p);

    // Readback and compare where interior lives
    let out_host = {
        let size = (n * std::mem::size_of::<f32>()) as u64;
        let read_buf = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read"),
            size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut enc = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("enc") });
        enc.copy_buffer_to_buffer(&out.buf, 0, &read_buf, 0, size);
        ctx.queue.submit(Some(enc.finish()));
        read_buf.slice(..).map_async(wgpu::MapMode::Read, |_| {});
        ctx.device.poll(wgpu::Maintain::Wait);
        let data = read_buf.slice(..).get_mapped_range();
        let mut v = vec![0.0f32; n];
        v.copy_from_slice(bytemuck::cast_slice(&data));
        drop(data);
        read_buf.unmap();
        v
    };

    // Compare in interior (skip halos)
    let mut max_abs = 0.0f32;
    for j in 0..height {
        for i in 0..width {
            let idx = (j + halo) * wt + (i + halo);
            let a_cpu = f_host[idx];
            let a_gpu = out_host[idx];
            max_abs = max_abs.max((a_cpu - a_gpu).abs());
        }
    }
    assert!(max_abs <= 1e-5, "max_abs = {}", max_abs);
}

#[test]
#[ignore]
fn single_v_cycle_reduces_residual() {
    // 128x128 interior + halo=2
    let width = 128usize;
    let height = 128usize;
    let halo = 2usize;
    let dx = 1.0f32;
    let d = 2.0f32;
    let k = 3.0f32;
    let (_w_true, f_host) = build_manufactured_rhs(width, height, halo, dx, d, k);
    let wt = width + 2 * halo;
    let ht = height + 2 * halo;
    let n = wt * ht;

    let ctx = pollster::block_on(engine::gpu::GpuContext::new());
    let mut flex = FlexGpu::new(&ctx);
    let dims = TileDims { width: width as u32, height: height as u32, halo: halo as u32, dx };
    let p = FlexParams { d, k, wj_omega: 0.7, nu1: 2, nu2: 2, levels: 1 };

    let mut w = GpuTex::new_storage(&ctx, n, "w");
    let f = {
        let t = GpuTex::new_storage(&ctx, n, "f");
        ctx.queue.write_buffer(&t.buf, 0, bytemuck::cast_slice(&f_host));
        t
    };
    // initialize w=0
    let zeros = vec![0.0f32; n];
    ctx.queue.write_buffer(&w.buf, 0, bytemuck::cast_slice(&zeros));
    let mut scratch = FlexScratch::new(&ctx, n);

    let stats = flex.v_cycle(&ctx, dims, &mut w, &f, &mut scratch, &p);
    assert!(
        stats.res_out < stats.res_in * 0.5,
        "residual not reduced enough: {} -> {}",
        stats.res_in,
        stats.res_out
    );
}
