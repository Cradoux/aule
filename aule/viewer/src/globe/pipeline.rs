use wgpu::util::DeviceExt;

use super::mesh::GlobeMesh;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Globals {
    view_proj: [[f32; 4]; 4],
    radius: f32,
    exagger: f32,
    debug_flags: u32,
    _pad: f32,
    d_max: f32,
    h_max: f32,
    height_scale: f32,
    _pad2: [f32; 2],
    _pad3: [f32; 3],
}

pub struct GlobeRenderer {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_lines: wgpu::RenderPipeline,
    pub bind_group0: wgpu::BindGroup,
    pub bind_group1: wgpu::BindGroup,
    pub uniform_buf: wgpu::Buffer,
    pub height_buf: wgpu::Buffer,
    pub lut_tex: wgpu::Texture,
    #[allow(dead_code)]
    pub lut_view: wgpu::TextureView,
    #[allow(dead_code)]
    pub lut_sampler: wgpu::Sampler,
}

impl GlobeRenderer {
    pub fn new(
        device: &wgpu::Device,
        surface_format: wgpu::TextureFormat,
        vertex_count: u32,
    ) -> Self {
        let globals_init = Globals {
            view_proj: [[0.0; 4]; 4],
            radius: 1.0,
            exagger: 0.0,
            debug_flags: 0,
            _pad: 0.0,
            d_max: 4000.0,
            h_max: 4000.0,
            height_scale: 1.0 / 6_371_000.0,
            _pad2: [0.0, 0.0],
            _pad3: [0.0, 0.0, 0.0],
        };
        let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("globe uniforms"),
            contents: bytemuck::bytes_of(&globals_init),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let height_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("globe height ssbo"),
            size: (vertex_count as usize * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // LUT texture identical to raster path (512x1 RGBA8UnormSrgb)
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

        let bgl0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("globe bgl0 uniforms"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        min_binding_size: None,
                        has_dynamic_offset: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D1,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });

        let bgl1 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("globe bgl1 height"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    min_binding_size: None,
                    has_dynamic_offset: false,
                },
                count: None,
            }],
        });

        let bind_group0 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("globe bg0"),
            layout: &bgl0,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: uniform_buf.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&lut_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::Sampler(&lut_sampler),
                },
            ],
        });

        let bind_group1 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("globe bg1"),
            layout: &bgl1,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: height_buf.as_entire_binding(),
            }],
        });

        let vert = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("globe vert"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../shaders/globe.vert.wgsl").into()),
        });

        let frag = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("globe frag"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../shaders/globe.frag.wgsl").into()),
        });

        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("globe pl"),
            bind_group_layouts: &[&bgl0, &bgl1],
            push_constant_ranges: &[],
        });

        let vertex_buffers = [wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<super::mesh::GlobeVertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                wgpu::VertexAttribute {
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                },
                wgpu::VertexAttribute {
                    shader_location: 1,
                    format: wgpu::VertexFormat::Uint32,
                    offset: 12,
                },
            ],
        }];

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("globe pipeline"),
            layout: Some(&pl),
            vertex: wgpu::VertexState {
                module: &vert,
                entry_point: "main",
                buffers: &vertex_buffers,
            },
            fragment: Some(wgpu::FragmentState {
                module: &frag,
                entry_point: "main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::LessEqual,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        // Simple line pipeline shares the same vertex layout and bg layouts
        let pl_lines = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("globe pl lines"),
            bind_group_layouts: &[&bgl0],
            push_constant_ranges: &[],
        });
        let shader_lines_src = r#"
struct Globals {
    view_proj : mat4x4<f32>,
    radius    : f32,
    exagger   : f32,
    debug_flags : u32,
    _pad : f32,
};
@group(0) @binding(0) var<uniform> G : Globals;

struct VSIn {
    @location(0) pos_unit : vec3<f32>,
    @location(1) vid      : u32,
};
struct VSOut { @builtin(position) pos_clip : vec4<f32> };

@vertex fn vmain(i: VSIn) -> VSOut {
    var o: VSOut;
    let world = normalize(i.pos_unit) * G.radius;
    o.pos_clip = G.view_proj * vec4<f32>(world, 1.0);
    return o;
}

@fragment fn fmain() -> @location(0) vec4<f32> {
    return vec4<f32>(1.0, 1.0, 1.0, 0.35);
}
"#;
        let shader_lines = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("globe lines shader"),
            source: wgpu::ShaderSource::Wgsl(shader_lines_src.into()),
        });
        let pipeline_lines = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("globe pipeline lines"),
            layout: Some(&pl_lines),
            vertex: wgpu::VertexState {
                module: &shader_lines,
                entry_point: "vmain",
                buffers: &vertex_buffers,
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader_lines,
                entry_point: "fmain",
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: false,
                depth_compare: wgpu::CompareFunction::LessEqual,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        Self {
            pipeline,
            pipeline_lines,
            bind_group0,
            bind_group1,
            uniform_buf,
            height_buf,
            lut_tex,
            lut_view,
            lut_sampler,
        }
    }

    pub fn upload_heights(&self, queue: &wgpu::Queue, heights: &[f32]) {
        // Render-only clamping for visualization; scientific arrays remain untouched
        let mut tmp: Vec<f32> = Vec::with_capacity(heights.len());
        for &v in heights {
            let z = if v.is_finite() { v } else { 0.0 };
            tmp.push(z.clamp(-11_000.0, 9_000.0));
        }
        queue.write_buffer(&self.height_buf, 0, bytemuck::cast_slice(&tmp));
    }

    pub fn update_uniforms(
        &self,
        queue: &wgpu::Queue,
        view_proj: [[f32; 4]; 4],
        radius: f32,
        exagger: f32,
        debug_flags: u32,
        d_max: f32,
        h_max: f32,
        height_scale: f32,
    ) {
        let u = Globals {
            view_proj,
            radius,
            exagger,
            debug_flags,
            _pad: 0.0,
            d_max,
            h_max,
            height_scale,
            _pad2: [0.0, 0.0],
            _pad3: [0.0, 0.0, 0.0],
        };
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    pub fn write_lut_from_overlay(&self, queue: &wgpu::Queue, ov: &crate::overlay::OverlayState) {
        let mut px: Vec<u8> = Vec::with_capacity(512 * 4);
        for i in 0..256 {
            let depth = (i as f32) / 255.0 * ov.hypso_d_max.max(1.0);
            let col = crate::overlay::ocean_color32(depth, ov.hypso_d_max.max(1.0));
            px.extend_from_slice(&[col.r(), col.g(), col.b(), 255]);
        }
        for i in 0..256 {
            let elev = (i as f32) / 255.0 * ov.hypso_h_max.max(1.0);
            let col =
                crate::overlay::land_color32(elev, ov.hypso_h_max.max(1.0), ov.hypso_snowline);
            px.extend_from_slice(&[col.r(), col.g(), col.b(), 255]);
        }
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

    pub fn draw<'a>(&'a self, rpass: &mut wgpu::RenderPass<'a>, mesh: &'a GlobeMesh) {
        rpass.set_pipeline(&self.pipeline);
        rpass.set_bind_group(0, &self.bind_group0, &[]);
        rpass.set_bind_group(1, &self.bind_group1, &[]);
        rpass.set_vertex_buffer(0, mesh.vertex_buf.slice(..));
        rpass.set_index_buffer(mesh.index_buf.slice(..), wgpu::IndexFormat::Uint32);
        rpass.draw_indexed(0..mesh.index_count, 0, 0..1);
    }

    pub fn draw_lines<'a>(&'a self, rpass: &mut wgpu::RenderPass<'a>, mesh: &'a GlobeMesh) {
        if mesh.line_buf.is_some() {
            rpass.set_pipeline(&self.pipeline_lines);
            rpass.set_bind_group(0, &self.bind_group0, &[]);
            rpass.set_vertex_buffer(0, mesh.vertex_buf.slice(..));
            if let Some(lb) = mesh.line_buf.as_ref() {
                rpass.set_index_buffer(lb.slice(..), wgpu::IndexFormat::Uint32);
            }
            rpass.draw_indexed(0..mesh.line_count, 0, 0..1);
        }
    }
}
