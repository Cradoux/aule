use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct GlobeVertex {
    pub pos_unit: [f32; 3],
    pub vid: u32,
}

pub struct GlobeMesh {
    pub vertex_buf: wgpu::Buffer,
    pub index_buf: wgpu::Buffer,
    pub index_count: u32,
    pub line_buf: Option<wgpu::Buffer>,
    pub line_count: u32,
    pub vertex_count: u32,
}

pub fn build_globe_mesh(device: &wgpu::Device, grid: &engine::grid::Grid) -> GlobeMesh {
    let vertex_count = grid.num_vertices() as u32;
    let mut verts: Vec<GlobeVertex> = Vec::with_capacity(vertex_count as usize);
    for (vid, p) in grid.pos_xyz.iter().enumerate() {
        verts.push(GlobeVertex { pos_unit: [p[0], p[1], p[2]], vid: vid as u32 });
    }

    // Indices: use per-face triangulation with shared vertex ids, exactly as Grid::new built tris
    let f = grid.freq() as usize;
    let mut indices: Vec<u32> = Vec::with_capacity((20 * f * f * 3) as usize);
    for face in 0..20usize {
        // DOWN triangles: [i,j],[i+1,j],[i,j+1]
        for i in 0..f {
            for j in 0..(f - i) {
                let a = grid.idx(face, i as u32, j as u32) as u32;
                let b = grid.idx(face, (i + 1) as u32, j as u32) as u32;
                let c = grid.idx(face, i as u32, (j + 1) as u32) as u32;
                indices.push(a);
                indices.push(b);
                indices.push(c);
            }
        }
        // UP triangles: [i+1,j],[i+1,j+1],[i,j+1]
        for i in 0..f {
            for j in 0..(f - i - 1) {
                let b = grid.idx(face, (i + 1) as u32, j as u32) as u32;
                let d = grid.idx(face, (i + 1) as u32, (j + 1) as u32) as u32;
                let c = grid.idx(face, i as u32, (j + 1) as u32) as u32;
                indices.push(b);
                indices.push(d);
                indices.push(c);
            }
        }
    }

    let vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("globe vertices"),
        contents: bytemuck::cast_slice(&verts),
        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
    });
    let index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("globe indices"),
        contents: bytemuck::cast_slice(&indices),
        usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
    });

    GlobeMesh {
        vertex_buf,
        index_buf,
        index_count: indices.len() as u32,
        line_buf: None,
        line_count: 0,
        vertex_count,
    }
}
