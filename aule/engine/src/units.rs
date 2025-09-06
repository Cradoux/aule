//! Minimal units-of-measure newtypes to make distances and times explicit.
//! Conversions are explicit; mixing units requires an intentional conversion.

/// Distance in meters (m).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Meters(pub f64);

/// Distance in kilometers (km).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Kilometers(pub f64);

/// Time in million years (Myr).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Myr(pub f64);

/// Velocity in meters per year (m/yr).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct MetersPerYear(pub f64);

impl Meters {
    /// Construct from a raw f64 in meters.
    pub fn new(v: f64) -> Self {
        Self(v)
    }
    /// Extract the raw numeric value in meters.
    pub fn value(self) -> f64 {
        self.0
    }
}

impl Kilometers {
    /// Construct from a raw f64 in kilometers.
    pub fn new(v: f64) -> Self {
        Self(v)
    }
    /// Extract the raw numeric value in kilometers.
    pub fn value(self) -> f64 {
        self.0
    }
}

impl Myr {
    /// Construct from a raw f64 in Myr.
    pub fn new(v: f64) -> Self {
        Self(v)
    }
    /// Extract the raw numeric value in Myr.
    pub fn value(self) -> f64 {
        self.0
    }
}

impl MetersPerYear {
    /// Construct from a raw f64 in m/yr.
    pub fn new(v: f64) -> Self {
        Self(v)
    }
    /// Extract the raw numeric value in m/yr.
    pub fn value(self) -> f64 {
        self.0
    }
}

// Explicit conversions
impl From<Kilometers> for Meters {
    fn from(km: Kilometers) -> Self {
        Meters(km.0 * 1000.0)
    }
}

// Helper conversion functions
/// Shorthand constructor for kilometers.
pub fn km(v: f64) -> Kilometers {
    Kilometers::new(v)
}
/// Shorthand constructor for meters.
pub fn m(v: f64) -> Meters {
    Meters::new(v)
}
/// Shorthand constructor for Myr.
pub fn myr(v: f64) -> Myr {
    Myr::new(v)
}
