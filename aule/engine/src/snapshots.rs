//! Snapshot writers for archiving simulation checkpoints (T-400).

use std::io::Write;

/// Write a CSV file containing the current depth field.
///
/// Format:
/// - First line: `# t_myr=<time>`
/// - Header: `index,depth_m`
/// - One row per cell with 0-based index and depth in meters (positive = down)
///
/// Errors are bubbled up from the filesystem.
pub fn write_csv_depth(path: &std::path::Path, t_myr: f64, depth_m: &[f32]) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "# t_myr={:.6}", t_myr)?;
    writeln!(file, "index,depth_m")?;
    for (i, &d) in depth_m.iter().enumerate() {
        // Clamp to finite numbers to avoid CSV pollution
        let v = if d.is_finite() { d } else { 0.0 };
        writeln!(file, "{},{}", i, v)?;
    }
    Ok(())
}
