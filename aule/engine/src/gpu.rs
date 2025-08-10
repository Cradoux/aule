//! Minimal GPU device/queue helper for tests (T-020).

use wgpu::{Device, Instance, Queue};

/// Minimal GPU context for tests and utilities.
pub struct GpuContext {
    /// Instance used to create adapters
    pub instance: Instance,
    /// Logical device
    pub device: Device,
    /// Submission queue
    pub queue: Queue,
}

impl GpuContext {
    /// Create a new GPU context using default instance and a high-performance adapter.
    pub async fn new() -> Self {
        let instance = wgpu::Instance::default();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .unwrap_or_else(|| panic!("no suitable GPU adapters"));
        // For tests we need storage buffers; request reasonable defaults.
        let required_limits = wgpu::Limits::default();
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("engine-device"),
                    required_features: wgpu::Features::empty(),
                    required_limits,
                },
                None,
            )
            .await
            .unwrap_or_else(|e| panic!("request device: {e}"));
        Self { instance, device, queue }
    }
}
