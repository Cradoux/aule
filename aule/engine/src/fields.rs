//! Field views and tile-aware accessors (index indirection only).
//! T-020: GPU field SoA buffers and trivial round-trip helpers.

use crate::grid::tile::TileView;
use bytemuck::{Pod, Zeroable};

/// Example of a structure-of-arrays storage for a scalar field.
#[derive(Debug)]
pub struct ScalarField {
    /// Global SoA storage (length equals number of grid cells)
    pub data: Vec<f32>,
}

impl ScalarField {
    /// Create a new scalar field with `len` cells initialized to 0.0
    pub fn new(len: usize) -> Self {
        Self { data: vec![0.0; len] }
    }

    /// Get a tile view as slices into the global storage by indices.
    /// Build a tile-aware view that carries index slices into the global storage.
    pub fn view<'a>(&'a self, tv: TileView<'a>) -> ScalarTileView<'a> {
        ScalarTileView { interior: tv.interior, halo: tv.halo, data: &self.data }
    }
}

/// A tile view for a scalar field that holds index slices and a reference to the data.
/// A borrowed tile view for a scalar field holding index slices and data reference.
pub struct ScalarTileView<'a> {
    /// Interior cell indices for this tile
    pub interior: &'a [u32],
    /// Halo cell indices (excludes interior)
    pub halo: &'a [u32],
    /// Reference to the global scalar data buffer
    pub data: &'a [f32],
}

// --------- T-020 GPU SoA buffers (draft) ---------

#[repr(C)]
#[derive(Copy, Clone, Debug, Zeroable, Pod)]
pub struct IdMaskU32(pub u32);

pub struct DeviceFields {
    pub cells: usize,
    // continuous SoA buffers
    pub h: wgpu::Buffer,           // f32
    pub th_c: wgpu::Buffer,        // f32
    pub age_ocean: wgpu::Buffer,   // f32
    pub s: wgpu::Buffer,           // f32 (salinity placeholder)
    pub v: wgpu::Buffer,           // f32 (velocity mag placeholder)
    pub plate_id: wgpu::Buffer,    // u32
    pub b: wgpu::Buffer,           // f32 (bathymetry/elevation)
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub bind_group: wgpu::BindGroup,
}

impl DeviceFields {
    pub fn new(device: &wgpu::Device, cells: usize) -> Self {
        fn buf(device: &wgpu::Device, size: usize, label: &str, usage: wgpu::BufferUsages) -> wgpu::Buffer {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size as u64,
                usage,
                mapped_at_creation: false,
            })
        }
        let bytes_f32 = cells * std::mem::size_of::<f32>();
        let bytes_u32 = cells * std::mem::size_of::<u32>();
        let usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let h = buf(device, bytes_f32, "h", usage);
        let th_c = buf(device, bytes_f32, "th_c", usage);
        let age_ocean = buf(device, bytes_f32, "age_ocean", usage);
        let s = buf(device, bytes_f32, "S", usage);
        let v = buf(device, bytes_f32, "V", usage);
        let plate_id = buf(device, bytes_u32, "plate_id", usage);
        let b = buf(device, bytes_f32, "B", usage);

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("fields-bgl"),
            entries: &[
                Self::entry(0), Self::entry(1), Self::entry(2), Self::entry(3), Self::entry(4), Self::entry(5), Self::entry(6),
            ],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("fields-bg"),
            layout: &bind_group_layout,
            entries: &[
                Self::bg_entry(0, &h),
                Self::bg_entry(1, &th_c),
                Self::bg_entry(2, &age_ocean),
                Self::bg_entry(3, &s),
                Self::bg_entry(4, &v),
                Self::bg_entry(5, &plate_id),
                Self::bg_entry(6, &b),
            ],
        });
        Self { cells, h, th_c, age_ocean, s, v, plate_id, b, bind_group_layout, bind_group }
    }

    fn entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
        wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        }
    }
    fn bg_entry<'a>(binding: u32, buffer: &'a wgpu::Buffer) -> wgpu::BindGroupEntry<'a> {
        wgpu::BindGroupEntry { binding, resource: buffer.as_entire_binding() }
    }

    pub fn resize(&mut self, device: &wgpu::Device, new_cells: usize) {
        if new_cells == self.cells { return; }
        *self = Self::new(device, new_cells);
    }
}
