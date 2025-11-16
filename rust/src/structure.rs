use std::collections::HashMap;

/// 3D coordinate vector
#[derive(Debug, Clone, Copy)]
pub struct Coordinate {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Coordinate {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Calculate Euclidean distance to another coordinate
    pub fn distance_to(&self, other: &Coordinate) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// Frame data: maps residue number to CA atom coordinate
pub type FrameData = HashMap<i32, Coordinate>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinate_distance() {
        let c1 = Coordinate::new(0.0, 0.0, 0.0);
        let c2 = Coordinate::new(3.0, 4.0, 0.0);
        assert_eq!(c1.distance_to(&c2), 5.0);
    }
}

