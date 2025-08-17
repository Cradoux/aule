#[repr(C)]
#[derive(Clone, Copy)]
pub struct Uniforms {
    pub width: u32,
    pub height: u32,
    pub f: u32,
    pub palette_mode: u32,
    pub _pad0: u32,
    pub d_max: f32,
    pub h_max: f32,
    pub snowline: f32,
    pub eta_m: f32,
}

pub struct RasterGpu {
    pub pipeline: wgpu::ComputePipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub uniforms: wgpu::Buffer,
    pub face_vert_ids: wgpu::Buffer,
    pub face_offsets: wgpu::Buffer,
    pub face_geom: wgpu::Buffer,
    pub vertex_values: wgpu::Buffer,
    #[allow(dead_code)]
    pub out_tex: wgpu::Texture,
    pub out_view: wgpu::TextureView,
    pub bind_group: wgpu::BindGroup,
    pub width: u32,
    pub height: u32,
    pub face_count: usize,
}

impl RasterGpu {
    pub fn new(
        device: &wgpu::Device,
        shader_src: &wgpu::ShaderModule,
        width: u32,
        height: u32,
    ) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("raster bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
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
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba8Unorm,
                        view_dimension: wgpu::TextureViewDimension::D2,
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("raster pl"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("raster pipeline"),
            layout: Some(&pipeline_layout),
            module: shader_src,
            entry_point: "main",
        });

        let out_tex = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("raster out tex"),
            size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8Unorm,
            usage: wgpu::TextureUsages::STORAGE_BINDING
                | wgpu::TextureUsages::TEXTURE_BINDING
                | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let out_view = out_tex.create_view(&wgpu::TextureViewDescriptor::default());

        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster uniforms"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let empty = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("raster bg"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: uniforms.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: wgpu::BindingResource::TextureView(&out_view),
                },
            ],
        });

        let empty2 = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty2"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let empty3 = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty3"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let empty4 = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty4"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Self {
            pipeline,
            bind_group_layout,
            uniforms,
            face_vert_ids: empty,
            face_offsets: empty2,
            face_geom: empty3,
            vertex_values: empty4,
            out_tex,
            out_view,
            bind_group,
            width,
            height,
            face_count: 0,
        }
    }

    fn pack_uniforms(u: &Uniforms) -> [u8; 40] {
        let mut bytes = [0u8; 40];
        bytes[0..4].copy_from_slice(&u.width.to_le_bytes());
        bytes[4..8].copy_from_slice(&u.height.to_le_bytes());
        bytes[8..12].copy_from_slice(&u.f.to_le_bytes());
        bytes[12..16].copy_from_slice(&u.palette_mode.to_le_bytes());
        bytes[16..20].copy_from_slice(&u._pad0.to_le_bytes());
        bytes[20..24].copy_from_slice(&u.d_max.to_le_bytes());
        bytes[24..28].copy_from_slice(&u.h_max.to_le_bytes());
        bytes[28..32].copy_from_slice(&u.snowline.to_le_bytes());
        bytes[32..36].copy_from_slice(&u.eta_m.to_le_bytes());
        // pad last 4 bytes to 16-byte alignment
        bytes[36..40].copy_from_slice(&0u32.to_le_bytes());
        bytes
    }

    pub fn upload_inputs(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        u: &Uniforms,
        face_vert_ids: Vec<u32>,
        face_offsets: Vec<u32>,
        face_geom: &[[f32; 4]],
        vertex_values: &[f32],
    ) {
        // Recreate buffers and upload via write_buffer (no extra deps)
        let u_bytes = Self::pack_uniforms(u);
        self.uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster uniforms upload"),
            size: u_bytes.len() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.uniforms, 0, &u_bytes);

        let mut ids_bytes: Vec<u8> = Vec::with_capacity(face_vert_ids.len() * 4);
        for v in &face_vert_ids {
            ids_bytes.extend_from_slice(&v.to_le_bytes());
        }
        self.face_vert_ids = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_vert_ids"),
            size: ids_bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.face_vert_ids, 0, &ids_bytes);

        let mut offs_bytes: Vec<u8> = Vec::with_capacity(face_offsets.len() * 4);
        for v in &face_offsets {
            offs_bytes.extend_from_slice(&v.to_le_bytes());
        }
        self.face_offsets = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_offsets"),
            size: offs_bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.face_offsets, 0, &offs_bytes);

        let mut geom_bytes: Vec<u8> = Vec::with_capacity(face_geom.len() * 16);
        for v in face_geom {
            for f in v {
                geom_bytes.extend_from_slice(&f.to_le_bytes());
            }
        }
        self.face_geom = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_geom"),
            size: geom_bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.face_geom, 0, &geom_bytes);

        let mut vals_bytes: Vec<u8> = Vec::with_capacity(vertex_values.len() * 4);
        for f in vertex_values {
            vals_bytes.extend_from_slice(&f.to_le_bytes());
        }
        self.vertex_values = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster vertex_values"),
            size: vals_bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.vertex_values, 0, &vals_bytes);
        self.face_count = face_offsets.len();
        self.bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("raster bg updated"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: self.uniforms.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.face_vert_ids.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.face_offsets.as_entire_binding(),
                },
                wgpu::BindGroupEntry { binding: 3, resource: self.face_geom.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: self.vertex_values.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: wgpu::BindingResource::TextureView(&self.out_view),
                },
            ],
        });
        let _ = queue; // not used after init here
    }

    pub fn dispatch(&self, device: &wgpu::Device, queue: &wgpu::Queue) {
        let gx = (self.width + 7) / 8;
        let gy = (self.height + 7) / 8;
        let mut enc = device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("raster enc") });
        {
            let mut cpass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("raster pass"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&self.pipeline);
            cpass.set_bind_group(0, &self.bind_group, &[]);
            cpass.dispatch_workgroups(gx, gy, 1);
        }
        queue.submit(std::iter::once(enc.finish()));
    }

    pub fn write_uniforms(&self, queue: &wgpu::Queue, u: &Uniforms) {
        let bytes = Self::pack_uniforms(u);
        queue.write_buffer(&self.uniforms, 0, &bytes);
    }

    pub fn write_vertex_values(&self, queue: &wgpu::Queue, vals: &[f32]) {
        let mut buf: Vec<u8> = Vec::with_capacity(vals.len() * 4);
        for f in vals {
            buf.extend_from_slice(&f.to_le_bytes());
        }
        queue.write_buffer(&self.vertex_values, 0, &buf);
    }
}
