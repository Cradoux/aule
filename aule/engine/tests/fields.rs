use engine::gpu::GpuContext;
use engine::fields::DeviceFields;

fn roundtrip_sizes() -> Vec<usize> { vec![1024, 1536, 2048] }

#[test]
fn device_fields_roundtrip_and_resize() {
    let ctx = pollster::block_on(GpuContext::new());
    for n in roundtrip_sizes() {
        let mut df = DeviceFields::new(&ctx.device, n);
        // Write patterns
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("enc") });
        let h_host: Vec<f32> = (0..n).map(|i| i as f32 + 0.5).collect();
        let pid_host: Vec<u32> = (0..n as u32).collect();
        // staging for upload
        let h_staging = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("h_staging"), contents: bytemuck::cast_slice(&h_host), usage: wgpu::BufferUsages::COPY_SRC,
        });
        let pid_staging = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("pid_staging"), contents: bytemuck::cast_slice(&pid_host), usage: wgpu::BufferUsages::COPY_SRC,
        });
        encoder.copy_buffer_to_buffer(&h_staging, 0, &df.h, 0, (n * std::mem::size_of::<f32>()) as u64);
        encoder.copy_buffer_to_buffer(&pid_staging, 0, &df.plate_id, 0, (n * std::mem::size_of::<u32>()) as u64);
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::Maintain::Wait);

        // Read back
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("rd") });
        let h_read = ctx.device.create_buffer(&wgpu::BufferDescriptor { label: Some("h_read"), size: (n * 4) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let pid_read = ctx.device.create_buffer(&wgpu::BufferDescriptor { label: Some("pid_read"), size: (n * 4) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        encoder.copy_buffer_to_buffer(&df.h, 0, &h_read, 0, (n * 4) as u64);
        encoder.copy_buffer_to_buffer(&df.plate_id, 0, &pid_read, 0, (n * 4) as u64);
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::Maintain::Wait);

        // Map and compare
        {
            let (s, r) = futures_intrusive::channel::shared::oneshot_channel();
            h_read.slice(..).map_async(wgpu::MapMode::Read, move |v| s.send(v).ok());
            pollster::block_on(r.receive());
            let data = h_read.slice(..).get_mapped_range();
            let out: &[f32] = bytemuck::cast_slice(&data);
            assert_eq!(out, &h_host[..]);
        }
        {
            let (s, r) = futures_intrusive::channel::shared::oneshot_channel();
            pid_read.slice(..).map_async(wgpu::MapMode::Read, move |v| s.send(v).ok());
            pollster::block_on(r.receive());
            let data = pid_read.slice(..).get_mapped_range();
            let out: &[u32] = bytemuck::cast_slice(&data);
            assert_eq!(out, &pid_host[..]);
        }

        // Resize grow then shrink
        df.resize(&ctx.device, n + 128);
        df.resize(&ctx.device, n);
    }
}


