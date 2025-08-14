//! Color atlas loaders and sampling (viewer-only).
//!
//! - `.pal` continuous palette with linear-sRGB interpolation
//! - `.cls` discrete classes for biome preview (lat/elev keyed)
//!
//! Assets are embedded with include_str! so there is no runtime IO.

/// A single color stop along a 1D value domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColorStop {
    /// Domain value (units depend on palette; meters for hypsometry)
    pub v: f32,
    /// sRGB 8-bit triplet
    pub rgb: [u8; 3],
}

/// Continuous palette: strictly increasing stops.
#[derive(Debug, Clone, PartialEq)]
pub struct Palette {
    pub stops: Vec<ColorStop>,
    pub vmin: f32,
    pub vmax: f32,
}

/// Discrete biome class keyed by latitude and elevation.
#[derive(Debug, Clone, PartialEq)]
pub struct BiomeClass {
    pub name: String,
    pub min_lat: f32,
    pub max_lat: f32,
    pub min_elev_m: f32,
    pub max_elev_m: f32,
    pub rgb: [u8; 3],
}

/// Embedded default assets
pub const HYPS_DEFAULT_STR: &str = include_str!("../assets/atlases/hyps_default.pal");
#[allow(dead_code)]
pub const BIOME_LAT_STR: &str = include_str!("../assets/atlases/biome_lat.cls");

/// Parse a `.pal` continuous palette.
///
/// Format: lines of `value  #RRGGBB`, blank lines and `#` comments are ignored.
pub fn parse_pal(src: &str) -> Result<Palette, String> {
    let mut stops: Vec<ColorStop> = Vec::new();
    for (lineno, raw) in src.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        // Remove inline comments after color (keep until end of color token)
        // Strategy: split by whitespace, find first token that starts with '#'
        let mut parts: Vec<&str> = Vec::new();
        let mut color_token: Option<&str> = None;
        for tok in line.split_whitespace() {
            if tok.starts_with('#') {
                color_token = Some(tok);
                break;
            } else {
                parts.push(tok);
            }
        }
        let v_str =
            parts.first().ok_or_else(|| format!(".pal: line {} missing value", lineno + 1))?;
        let v: f32 = v_str
            .parse::<f32>()
            .map_err(|_| format!(".pal: line {} bad value '{}", lineno + 1, v_str))?;
        let col_tok =
            color_token.ok_or_else(|| format!(".pal: line {} missing color", lineno + 1))?;
        let rgb = parse_hex_rgb(col_tok).map_err(|e| format!(".pal: line {} {}", lineno + 1, e))?;
        stops.push(ColorStop { v, rgb });
    }
    if stops.is_empty() {
        return Err(".pal: no stops".to_string());
    }
    // Ensure strictly increasing v
    for i in 1..stops.len() {
        if stops[i].v.partial_cmp(&stops[i - 1].v) != Some(std::cmp::Ordering::Greater) {
            return Err(format!(".pal: values must be strictly increasing at index {}", i));
        }
    }
    let vmin = stops[0].v;
    let vmax = stops[stops.len() - 1].v;
    Ok(Palette { stops, vmin, vmax })
}

/// Parse a `.cls` discrete classes file.
///
/// Format: rows of `name  min_lat  max_lat  min_elev_m  max_elev_m  #RRGGBB`
#[allow(dead_code)]
pub fn parse_cls(src: &str) -> Result<Vec<BiomeClass>, String> {
    let mut classes: Vec<BiomeClass> = Vec::new();
    for (lineno, raw) in src.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        // Split by whitespace; expect 6 tokens
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 6 {
            return Err(format!(
                ".cls: line {} expected 6 columns, got {}",
                lineno + 1,
                toks.len()
            ));
        }
        let name = toks[0].to_string();
        let min_lat: f32 =
            toks[1].parse().map_err(|_| format!(".cls: line {} bad min_lat", lineno + 1))?;
        let max_lat: f32 =
            toks[2].parse().map_err(|_| format!(".cls: line {} bad max_lat", lineno + 1))?;
        let min_elev_m: f32 =
            toks[3].parse().map_err(|_| format!(".cls: line {} bad min_elev_m", lineno + 1))?;
        let max_elev_m: f32 =
            toks[4].parse().map_err(|_| format!(".cls: line {} bad max_elev_m", lineno + 1))?;
        let rgb = parse_hex_rgb(toks[5]).map_err(|e| format!(".cls: line {} {}", lineno + 1, e))?;
        classes.push(BiomeClass { name, min_lat, max_lat, min_elev_m, max_elev_m, rgb });
    }
    Ok(classes)
}

#[inline]
fn parse_hex_rgb(tok: &str) -> Result<[u8; 3], String> {
    let s = tok.trim();
    if !s.starts_with('#') {
        return Err("expected #RRGGBB".to_string());
    }
    let hex = &s[1..];
    if hex.len() != 6 {
        return Err("expected 6 hex digits".to_string());
    }
    let r = u8::from_str_radix(&hex[0..2], 16).map_err(|_| "bad R".to_string())?;
    let g = u8::from_str_radix(&hex[2..4], 16).map_err(|_| "bad G".to_string())?;
    let b = u8::from_str_radix(&hex[4..6], 16).map_err(|_| "bad B".to_string())?;
    Ok([r, g, b])
}

#[inline]
#[allow(dead_code)]
fn srgb_u8_to_linear(rgb: [u8; 3]) -> [f32; 3] {
    [srgb_to_linear(rgb[0]), srgb_to_linear(rgb[1]), srgb_to_linear(rgb[2])]
}

#[inline]
#[allow(dead_code)]
fn srgb_to_linear(c: u8) -> f32 {
    let x = (c as f32) / 255.0;
    if x <= 0.04045 {
        x / 12.92
    } else {
        ((x + 0.055) / 1.055).powf(2.4)
    }
}

#[inline]
#[allow(dead_code)]
fn linear_to_srgb_u8(c: f32) -> u8 {
    let y = if c <= 0.003_130_8 { 12.92 * c } else { 1.055 * c.powf(1.0 / 2.4) - 0.055 };
    (y.clamp(0.0, 1.0) * 255.0 + 0.5).floor() as u8
}

/// Sample palette at x with linear-RGB interpolation (gamma-correct).
/// Values are clamped to [vmin, vmax].
#[allow(dead_code)]
pub fn sample_linear_srgb(p: &Palette, x: f32) -> [u8; 3] {
    let n = p.stops.len();
    if n == 0 {
        return [0, 0, 0];
    }
    if n == 1 {
        return p.stops[0].rgb;
    }
    let x = x.clamp(p.vmin, p.vmax);
    // Find segment [i, i+1] such that v_i <= x <= v_{i+1}
    let mut i = 0usize;
    if x >= p.stops[n - 1].v {
        i = n - 2;
    } else {
        // linear scan is fine for small N
        for k in 0..n - 1 {
            if x <= p.stops[k + 1].v {
                i = k;
                break;
            }
        }
    }
    let a = p.stops[i];
    let b = p.stops[i + 1];
    let t = if b.v > a.v { (x - a.v) / (b.v - a.v) } else { 0.0 };
    let la = srgb_u8_to_linear(a.rgb);
    let lb = srgb_u8_to_linear(b.rgb);
    let lr =
        [la[0] + t * (lb[0] - la[0]), la[1] + t * (lb[1] - la[1]), la[2] + t * (lb[2] - la[2])];
    [linear_to_srgb_u8(lr[0]), linear_to_srgb_u8(lr[1]), linear_to_srgb_u8(lr[2])]
}

use std::sync::OnceLock;

#[allow(dead_code)]
static HYPS_DEFAULT: OnceLock<Palette> = OnceLock::new();
#[allow(dead_code)]
static BIOME_DEFAULT: OnceLock<Vec<BiomeClass>> = OnceLock::new();

#[allow(dead_code)]
pub fn hyps_default_palette() -> &'static Palette {
    HYPS_DEFAULT.get_or_init(|| {
        parse_pal(HYPS_DEFAULT_STR).unwrap_or_else(|e| panic!("hyps_default.pal: {e}"))
    })
}

#[allow(dead_code)]
pub fn biome_default_classes() -> &'static [BiomeClass] {
    BIOME_DEFAULT
        .get_or_init(|| parse_cls(BIOME_LAT_STR).unwrap_or_else(|e| panic!("biome_lat.cls: {e}")))
}

/// Pick the first matching class by latitude and elevation.
#[allow(dead_code)]
pub fn pick_biome_color(lat_rad: f32, elev_m: f32) -> [u8; 3] {
    let lat_deg = lat_rad.to_degrees();
    for cls in biome_default_classes() {
        if lat_deg >= cls.min_lat
            && lat_deg <= cls.max_lat
            && elev_m >= cls.min_elev_m
            && elev_m <= cls.max_elev_m
        {
            return cls.rgb;
        }
    }
    // Fallback: gray
    [128, 128, 128]
}
