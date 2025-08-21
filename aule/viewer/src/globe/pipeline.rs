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
}

pub struct GlobeRenderer {
    pub pipeline: wgpu::RenderPipeline,
    pub bind_group0: wgpu::BindGroup,
    pub bind_group1: wgpu::BindGroup,
    pub uniform_buf: wgpu::Buffer,
    pub height_buf: wgpu::Buffer,
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

        let bgl0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("globe bgl0 uniforms"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    min_binding_size: None,
                    has_dynamic_offset: false,
                },
                count: None,
            }],
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
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buf.as_entire_binding(),
            }],
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
                ..Default::default()
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        Self { pipeline, bind_group0, bind_group1, uniform_buf, height_buf }
    }

    pub fn upload_heights(&self, queue: &wgpu::Queue, heights: &[f32]) {
        queue.write_buffer(&self.height_buf, 0, bytemuck::cast_slice(heights));
    }

    pub fn update_uniforms(
        &self,
        queue: &wgpu::Queue,
        view_proj: [[f32; 4]; 4],
        radius: f32,
        exagger: f32,
        debug_flags: u32,
    ) {
        let u = Globals { view_proj, radius, exagger, debug_flags, _pad: 0.0 };
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    pub fn draw<'a>(&'a self, rpass: &mut wgpu::RenderPass<'a>, mesh: &'a GlobeMesh) {
        rpass.set_pipeline(&self.pipeline);
        rpass.set_bind_group(0, &self.bind_group0, &[]);
        rpass.set_bind_group(1, &self.bind_group1, &[]);
        rpass.set_vertex_buffer(0, mesh.vertex_buf.slice(..));
        rpass.set_index_buffer(mesh.index_buf.slice(..), wgpu::IndexFormat::Uint32);
        rpass.draw_indexed(0..mesh.index_count, 0, 0..1);
    }
}
