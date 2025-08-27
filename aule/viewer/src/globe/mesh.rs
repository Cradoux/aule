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
    let f = grid.freq();
    let fz = f as usize;
    let mut indices: Vec<u32> = Vec::with_capacity((20 * fz * fz * 3) as usize);
    let (face_vert_ids, face_offsets) = grid.face_vertex_table();
    let row_base = |ii: u32| -> u32 { ii * (f + 1) - (ii * ii.saturating_sub(1)) / 2 };
    let idx_at = |face: usize, i: u32, j: u32| -> u32 {
        let base = face_offsets[face] + row_base(i) + j;
        face_vert_ids[base as usize]
    };
    for face in 0..20usize {
        // DOWN triangles: [i,j],[i+1,j],[i,j+1]
        for i in 0..fz {
            for j in 0..(fz - i) {
                let a = idx_at(face, i as u32, j as u32);
                let b = idx_at(face, (i + 1) as u32, j as u32);
                let c = idx_at(face, i as u32, (j + 1) as u32);
                indices.push(a);
                indices.push(b);
                indices.push(c);
            }
        }
        // UP triangles: [i+1,j],[i+1,j+1],[i,j+1]
        for i in 0..fz {
            for j in 0..(fz - i - 1) {
                let b = idx_at(face, (i + 1) as u32, j as u32);
                let d = idx_at(face, (i + 1) as u32, (j + 1) as u32);
                let c = idx_at(face, i as u32, (j + 1) as u32);
                indices.push(b);
                indices.push(d);
                indices.push(c);
            }
        }
    }

    // Build unique edge set for optional wireframe
    use std::collections::HashSet;
    let mut edges: HashSet<(u32, u32)> = HashSet::new();
    // Iterate tris again using the same local accessor
    for face in 0..20usize {
        for i in 0..fz {
            for j in 0..(fz - i) {
                let a = idx_at(face, i as u32, j as u32);
                let b = idx_at(face, (i + 1) as u32, j as u32);
                let c = idx_at(face, i as u32, (j + 1) as u32);
                let mut add = |u: u32, v: u32| {
                    let e = if u < v { (u, v) } else { (v, u) };
                    edges.insert(e);
                };
                add(a, b);
                add(b, c);
                add(c, a);
            }
        }
        for i in 0..fz {
            for j in 0..(fz - i - 1) {
                let b = idx_at(face, (i + 1) as u32, j as u32);
                let d = idx_at(face, (i + 1) as u32, (j + 1) as u32);
                let c = idx_at(face, i as u32, (j + 1) as u32);
                let mut add = |u: u32, v: u32| {
                    let e = if u < v { (u, v) } else { (v, u) };
                    edges.insert(e);
                };
                add(b, d);
                add(d, c);
                add(c, b);
            }
        }
    }
    let mut line_list: Vec<u32> = Vec::with_capacity(edges.len() * 2);
    for (u, v) in edges.into_iter() {
        line_list.push(u);
        line_list.push(v);
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
    let line_buf = if !line_list.is_empty() {
        Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("globe lines"),
            contents: bytemuck::cast_slice(&line_list),
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
        }))
    } else {
        None
    };

    GlobeMesh {
        vertex_buf,
        index_buf,
        index_count: indices.len() as u32,
        line_buf,
        line_count: line_list.len() as u32,
        vertex_count,
    }
}
