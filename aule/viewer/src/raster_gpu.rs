#[repr(C)]
#[derive(Clone, Copy)]
pub struct Uniforms {
    pub width: u32,
    pub height: u32,
    pub f: u32,
    pub palette_mode: u32,
    pub debug_flags: u32,
    pub d_max: f32,
    pub h_max: f32,
    pub snowline: f32,
    pub eta_m: f32,
    pub inv_dmax: f32,
    pub inv_hmax: f32,
}

pub struct RasterGpu {
    pub pipeline: wgpu::ComputePipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub uniforms: wgpu::Buffer,
    pub face_vert_ids: wgpu::Buffer,
    pub face_offsets: wgpu::Buffer,
    pub face_geom: wgpu::Buffer,
    pub vertex_values: wgpu::Buffer,
    pub face_tri_offsets: wgpu::Buffer,
    pub face_tri_indices: wgpu::Buffer,
    pub face_edge_info: wgpu::Buffer,
    pub debug_face_tri: wgpu::Buffer,
    pub cpu_face_pick: wgpu::Buffer,
    pub lut_tex: wgpu::Texture,
    pub lut_view: wgpu::TextureView,
    pub lut_sampler: wgpu::Sampler,
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
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D1,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 6,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 7,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba8Unorm,
                        view_dimension: wgpu::TextureViewDimension::D2,
                    },
                    count: None,
                },
                // Drop bindings 8,9 to stay under storage buffer limits
                wgpu::BindGroupLayoutEntry {
                    binding: 10,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 11,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 12,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
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

        // LUT texture (1D x 512)
        let lut_tex = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("palette LUT"),
            size: wgpu::Extent3d { width: 512, height: 1, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D1,
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
            view_formats: &[],
        });
        let lut_view = lut_tex.create_view(&wgpu::TextureViewDescriptor {
            label: Some("palette LUT view"),
            format: Some(wgpu::TextureFormat::Rgba8UnormSrgb),
            dimension: Some(wgpu::TextureViewDimension::D1),
            aspect: wgpu::TextureAspect::All,
            base_mip_level: 0,
            mip_level_count: None,
            base_array_layer: 0,
            array_layer_count: None,
        });
        let lut_sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("palette LUT sampler"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::FilterMode::Nearest,
            ..Default::default()
        });

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
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&lut_view),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: wgpu::BindingResource::Sampler(&lut_sampler),
                },
                wgpu::BindGroupEntry { binding: 5, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: wgpu::BindingResource::TextureView(&out_view),
                },
                // bindings 8 and 9 unused in the layout; use placeholders here
                wgpu::BindGroupEntry { binding: 10, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 11, resource: empty.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 12, resource: empty.as_entire_binding() },
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
        let empty5 = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty5"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let empty6 = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("empty6"),
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
            face_tri_offsets: empty5,
            face_tri_indices: empty6,
            face_edge_info: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("empty7"),
                size: 4,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            debug_face_tri: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("empty8"),
                size: 4,
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }),
            cpu_face_pick: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("empty9"),
                size: 4,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            lut_tex,
            lut_view,
            lut_sampler,
            out_tex,
            out_view,
            bind_group,
            width,
            height,
            face_count: 0,
        }
    }

    fn pack_uniforms(u: &Uniforms) -> [u8; 56] {
        let mut bytes = [0u8; 56];
        bytes[0..4].copy_from_slice(&u.width.to_le_bytes());
        bytes[4..8].copy_from_slice(&u.height.to_le_bytes());
        bytes[8..12].copy_from_slice(&u.f.to_le_bytes());
        bytes[12..16].copy_from_slice(&u.palette_mode.to_le_bytes());
        bytes[16..20].copy_from_slice(&u.debug_flags.to_le_bytes());
        bytes[20..24].copy_from_slice(&u.d_max.to_le_bytes());
        bytes[24..28].copy_from_slice(&u.h_max.to_le_bytes());
        bytes[28..32].copy_from_slice(&u.snowline.to_le_bytes());
        bytes[32..36].copy_from_slice(&u.eta_m.to_le_bytes());
        bytes[36..40].copy_from_slice(&u.inv_dmax.to_le_bytes());
        bytes[40..44].copy_from_slice(&u.inv_hmax.to_le_bytes());
        // pad to 16-byte alignment
        bytes[44..48].copy_from_slice(&0u32.to_le_bytes());
        bytes[48..52].copy_from_slice(&0u32.to_le_bytes());
        bytes[52..56].copy_from_slice(&0u32.to_le_bytes());
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
        face_edge_info: &[u32],
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

        // Build per-face triangle indices and offsets (lower then upper)
        let f = u.f;
        let faces_n = 20u32;
        let mut tri_offs: Vec<u32> = Vec::with_capacity(faces_n as usize);
        let mut tris: Vec<u32> = Vec::new();
        // helper to get global vertex id for (face,i,j)
        let row_base = |ii: u32| -> u32 {
            if ii == 0 {
                0
            } else {
                ii * (f + 1) - (ii * (ii - 1)) / 2
            }
        };
        for face in 0..faces_n {
            let f_off = face_offsets[face as usize];
            // triangular row linearization: rows iv = 0..F-1, iu = 0..(F-1-iv)
            tri_offs.push((tris.len() as u32) / 3);
            for iv in 0..f {
                let max_u = f - 1 - iv;
                for iu in 0..=max_u {
                    // lower triangle (iu,iv) → (iu+1,iv) → (iu,iv+1)
                    let id_ij = face_vert_ids[(f_off + row_base(iu) + iv) as usize];
                    let id_i1j = face_vert_ids[(f_off + row_base(iu + 1) + iv) as usize];
                    let id_ij1 = face_vert_ids[(f_off + row_base(iu) + (iv + 1)) as usize];
                    tris.push(id_ij);
                    tris.push(id_i1j);
                    tris.push(id_ij1);
                    // upper triangle exists when iu < max_u
                    if iu < max_u {
                        let id_i1j1 = face_vert_ids[(f_off + row_base(iu + 1) + (iv + 1)) as usize];
                        tris.push(id_i1j);
                        tris.push(id_i1j1);
                        tris.push(id_ij1);
                    }
                }
            }
        }
        // Upload tri buffers
        let mut tri_offs_bytes: Vec<u8> = Vec::with_capacity(tri_offs.len() * 4);
        for v in &tri_offs {
            tri_offs_bytes.extend_from_slice(&v.to_le_bytes());
        }
        // Keep buffers allocated but unused to maintain struct fields
        self.face_tri_offsets = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_tri_offsets (unused)"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.face_tri_indices = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_tri_indices (unused)"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Upload neighbor edge info (20 faces x 3 edges x 3 u32)
        let mut e_bytes: Vec<u8> = Vec::with_capacity(face_edge_info.len() * 4);
        for v in face_edge_info {
            e_bytes.extend_from_slice(&v.to_le_bytes());
        }
        self.face_edge_info = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster face_edge_info"),
            size: e_bytes.len() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&self.face_edge_info, 0, &e_bytes);

        // Allocate debug face/tri buffer: 2 u32 per pixel
        let dbg_len = (self.width as usize) * (self.height as usize) * 2 * 4;
        self.debug_face_tri = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster debug_face_tri"),
            size: dbg_len as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        // Allocate cpu_face_pick buffer: 1 u32 per pixel
        let cpu_len = (self.width as usize) * (self.height as usize) * 4;
        self.cpu_face_pick = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster cpu_face_pick"),
            size: cpu_len as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

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
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&self.lut_view),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: wgpu::BindingResource::Sampler(&self.lut_sampler),
                },
                wgpu::BindGroupEntry { binding: 5, resource: self.face_geom.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: self.vertex_values.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: wgpu::BindingResource::TextureView(&self.out_view),
                },
                // bindings 8,9 unused
                wgpu::BindGroupEntry {
                    binding: 10,
                    resource: self.face_edge_info.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 11,
                    resource: self.debug_face_tri.as_entire_binding(),
                },
                // binding 12 only when using CPU face pick; otherwise bind a small dummy to keep layout consistent
                wgpu::BindGroupEntry {
                    binding: 12,
                    resource: self.cpu_face_pick.as_entire_binding(),
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

    pub fn write_lut_from_overlay(&self, queue: &wgpu::Queue, ov: &crate::overlay::OverlayState) {
        // Build 512 RGBA rows: 0..255 ocean (deep->0), 256..511 land
        let mut px: Vec<u8> = Vec::with_capacity(512 * 4);
        for i in 0..256 {
            let depth = (i as f32) / 255.0 * ov.hypso_d_max.max(1.0);
            let col = crate::overlay::ocean_color32(depth, ov.hypso_d_max.max(1.0));
            px.push(col.r());
            px.push(col.g());
            px.push(col.b());
            px.push(255);
        }
        for i in 0..256 {
            let elev = (i as f32) / 255.0 * ov.hypso_h_max.max(1.0);
            let col =
                crate::overlay::land_color32(elev, ov.hypso_h_max.max(1.0), ov.hypso_snowline);
            px.push(col.r());
            px.push(col.g());
            px.push(col.b());
            px.push(255);
        }
        // NonZeroU32::new returns Some for non-zero; since our inputs are compile-time non-zero, bypass fallbacks
        let bpr = std::num::NonZeroU32::new(512 * 4).map(|nz| nz.into()).unwrap_or(512 * 4);
        let rpi = std::num::NonZeroU32::new(1).map(|nz| nz.into()).unwrap_or(1);
        let layout = wgpu::ImageDataLayout {
            offset: 0,
            bytes_per_row: Some(bpr),
            rows_per_image: Some(rpi),
        };
        queue.write_texture(
            wgpu::ImageCopyTexture {
                texture: &self.lut_tex,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            &px,
            layout,
            wgpu::Extent3d { width: 512, height: 1, depth_or_array_layers: 1 },
        );
    }

    pub fn write_cpu_face_pick(&self, queue: &wgpu::Queue, data: &[u32]) {
        let mut bytes: Vec<u8> = Vec::with_capacity(data.len() * 4);
        for v in data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        queue.write_buffer(&self.cpu_face_pick, 0, &bytes);
    }

    /// Read back the debug face/tri buffer written by the compute shader when debug bit 6 is set.
    /// Returns a Vec<u32> with length width*height*2 containing [face, tri] per pixel in row-major order.
    pub fn read_debug_face_tri(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
    ) -> Option<Vec<u32>> {
        let bytes_len = (self.width as usize) * (self.height as usize) * 2 * 4;
        if bytes_len == 0 {
            return None;
        }
        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("raster debug readback"),
            size: bytes_len as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("copy debug buffer"),
        });
        enc.copy_buffer_to_buffer(&self.debug_face_tri, 0, &staging, 0, bytes_len as u64);
        queue.submit(std::iter::once(enc.finish()));

        // Map and wait
        let slice = staging.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = sender.send(r.is_ok());
        });
        device.poll(wgpu::Maintain::Wait);
        if receiver.recv().ok() != Some(true) {
            return None;
        }
        let data = slice.get_mapped_range();
        let mut out: Vec<u32> = vec![0; bytes_len / 4];
        // Copy bytes into u32 vec (little endian)
        let mut i = 0usize;
        for chunk in data.chunks_exact(4) {
            out[i] = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            i += 1;
        }
        drop(data);
        staging.unmap();
        Some(out)
    }
}
