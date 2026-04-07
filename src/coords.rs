use crate::RADIUS_EARTH;
use std::f64::consts::PI;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Div, Mul, Sub};

/// Spherical coordinate system (radius, colatitude, longitude).
///
/// Angles are in radians.
#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct Spherical {
    /// Radius
    pub(crate) r: f64,

    /// Co-latitude in radians
    pub(crate) t: f64,

    /// Longitude in radians
    pub(crate) p: f64,
}
impl Spherical {
    pub(crate) fn to_cartesian(self) -> Cartesian {
        let sq = self.r * self.t.sin();
        let x = sq * self.p.cos();
        let y = sq * self.p.sin();
        let z = self.r * self.t.cos();

        Cartesian { x, y, z }
    }

    /// Latitude in degrees
    pub(crate) fn lat(&self) -> f64 {
        90.0 - self.t.to_degrees()
    }

    /// Longitude in degrees
    pub(crate) fn lon(&self) -> f64 {
        self.p.to_degrees()
    }

    /// Converts a magnetic field vector at `self` into Cartesian components.
    ///
    /// Coordinates of `b_field` are the radial, southward, and eastward components of the
    /// magnetic field, in units of nT.
    pub(crate) fn vec_to_cartesian(&self, b_field: SphericalVec) -> Cartesian {
        let st = self.t.sin();
        let ct = self.t.cos();
        let sp = self.p.sin();
        let cp = self.p.cos();
        let be = b_field.radial() * st + b_field.colatitude() * ct;

        Cartesian {
            x: be * cp - b_field.longitude() * sp,
            y: be * sp + b_field.longitude() * cp,
            z: b_field.radial() * ct - b_field.colatitude() * st,
        }
    }
}
impl Display for Spherical {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "(r: {}, t: {}, p: {})",
            self.r,
            90. - self.t.to_degrees(),
            self.p.to_degrees()
        )
    }
}

/// Cartesian coordinate system. Units of km.
#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct Cartesian {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

impl Cartesian {
    /// Convert to ['Spherical'] coordinates.
    ///
    /// **Note**: at the poles (x=0 and y=0) it is assumed that phi=0
    pub(crate) fn to_spherical(self) -> Spherical {
        let mut sq = self.x * self.x + self.y * self.y;
        let radius = (sq + self.z * self.z).sqrt();
        let theta: f64;
        let mut phi: f64;

        if sq == 0. {
            phi = 0.;
            theta = if self.z < 0.0 { PI } else { 0. };
        } else {
            sq = sq.sqrt();
            phi = self.y.atan2(self.x);
            theta = sq.atan2(self.z);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }
        }

        Spherical {
            r: radius,
            t: theta,
            p: phi,
        }
    }

    /// Calculates the magnitude of the vector.
    pub(crate) fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub(crate) fn to_geo(self) -> Geographic {
        Geographic { coords: self }
    }
}
impl Add for Cartesian {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
impl Sub for Cartesian {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}
impl Mul<f64> for Cartesian {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}
impl Mul<Cartesian> for f64 {
    type Output = Cartesian;

    fn mul(self, rhs: Cartesian) -> Self::Output {
        Self::Output {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl Div<f64> for Cartesian {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

/// Magnetic coordinates, derived from IGRF models.
#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct Magnetic {
    pub(crate) coords: Cartesian,
}

impl Magnetic {
    pub(crate) fn to_spherical(self) -> Spherical {
        self.coords.to_spherical()
    }
}
impl Display for Magnetic {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.coords.to_spherical())
    }
}

/// Vector represented by radial, colatitudinal, and longitudinal components.
#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct SphericalVec {
    pub(crate) coords: Cartesian,
}
impl SphericalVec {
    pub(crate) fn radial(&self) -> f64 {
        self.coords.x
    }
    pub(crate) fn colatitude(&self) -> f64 {
        self.coords.y
    }
    pub(crate) fn longitude(&self) -> f64 {
        self.coords.z
    }
}

/// Geographic coordinates, in a Cartesian basis.
#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct Geographic {
    pub(crate) coords: Cartesian,
}

impl Geographic {
    /// Calculates the magnitude of the vector.
    pub(crate) fn norm(&self) -> f64 {
        self.coords.norm()
    }
}
impl Display for Geographic {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.coords.to_spherical())
    }
}

#[derive(Debug)]
pub(crate) struct Geodetic {
    /// Radial component represents altitude above sea level, in km.
    pub(crate) coords: Spherical,
}
impl Geodetic {
    pub(crate) fn to_geocentric(&self) -> Geocentric {
        if self.coords.r < 10.0 || self.coords.r > 3000.0 {
            panic!("Radius not in km above sea level")
        }

        let a = 6378.1370; /* semi-major axis */
        let f = 1. / 298.257223563; /* flattening */
        let b = a * (1.0 - f); /* semi-minor axis */
        let a2 = a * a;
        let b2 = b * b;
        let theta = self.coords.t; /* colatitude in radians   */
        let st = theta.sin();
        let ct = theta.cos();
        let one = a2 * st * st;
        let two = b2 * ct * ct;
        let three = one + two;
        let rho = three.sqrt(); /* [km] */
        let r =
            (self.coords.r * (self.coords.r + 2.0 * rho) + (a2 * one + b2 * two) / three).sqrt(); /* [km] */
        let cd = (self.coords.r + rho) / r;
        let sd = (a2 - b2) / rho * ct * st / r;

        let coords = Spherical {
            r: r / RADIUS_EARTH,
            t: (ct * cd - st * sd).acos(),
            p: self.coords.p,
        };
        Geocentric { coords }
    }
}

/// Geocentric coordinate system.
///
/// Radial component is given as distance from the center of the Earth, in units of Earth's radius.
/// Angular components are
pub(crate) struct Geocentric {
    /// Radial component represents distance from center of the Earth, in units of Earth's radius.
    pub(crate) coords: Spherical,
}

impl Geocentric {
    /// Latitude in degrees
    pub(crate) fn lat(&self) -> f64 {
        self.coords.lat()
    }

    /// Longitude in degrees
    pub(crate) fn lon(&self) -> f64 {
        self.coords.lon()
    }

    /// Altitude above sea level (assuming spherical Earth)
    pub(crate) fn alt(&self) -> f64 {
        (self.coords.r - 1.0) * RADIUS_EARTH
    }

    pub(crate) fn to_geodetic(&self) -> Geodetic {
        let a = 6378.1370; /* semi-major axis */
        let f = 1.0 / 298.257223563; /* flattening */
        let ee = (2.0 - f) * f;
        let e4 = ee * ee;
        let aa = a * a;

        let theta = self.coords.t;
        let phi = self.coords.p;

        let st = theta.sin();
        let ct = theta.cos();
        let sp = phi.sin();
        let cp = phi.cos();

        let x = self.coords.r * RADIUS_EARTH * st * cp;
        let y = self.coords.r * RADIUS_EARTH * st * sp;
        let z = self.coords.r * RADIUS_EARTH * ct;

        let k0i = 1.0 - ee;
        let pp = x * x + y * y;
        let zeta = k0i * z * z / aa;
        let rho = (pp / aa + zeta - e4) / 6.0;
        let s = e4 * zeta * pp / (4.0 * aa);
        let rho3 = rho * rho * rho;
        let t = (rho3 + s + (s * (s + 2.0 * rho3)).sqrt()).powf(1.0 / 3.0);
        let u = rho + t + rho * rho / t;
        let v = (u * u + e4 * zeta).sqrt();
        let w = ee * (u + v - zeta) / (2.0 * v);
        let kappa = 1.0 + ee * ((u + v + w * w).sqrt() + w) / (u + v);

        let coords = Spherical {
            r: (pp + z * z * kappa * kappa).sqrt() / ee * (1.0 / kappa - k0i),
            t: (z * kappa).atan2(pp.sqrt()),
            p: phi,
        };
        Geodetic { coords }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_car2sph() {
        let xyz = Cartesian {
            x: 0.759081899,
            y: -0.330058202,
            z: 0.837520183,
        };
        let rtp = xyz.to_spherical();
        assert_relative_eq!(rtp.r, 1.177532931, max_relative = 1.0e-8);
        assert_relative_eq!(rtp.t, 0.779521604, max_relative = 1.0e-8);
        assert_relative_eq!(rtp.p, 5.873032935, max_relative = 1.0e-8);
    }

    #[test]
    fn test_spherical() {
        let rtp = Spherical {
            r: 1.177532882,
            t: 0.779521585,
            p: -0.410152374,
        };
        let xyz = rtp.to_cartesian();
        assert_relative_eq!(xyz.x, 0.759081852, max_relative = 1.0e-8);
        assert_relative_eq!(xyz.y, -0.330058183, max_relative = 1.0e-8);
        assert_relative_eq!(xyz.z, 0.837520163, max_relative = 1.0e-8);

        let rtp = Spherical {
            r: 1.0,
            t: 0.779521583,
            p: 5.873032931,
        };
        let b_vec = Cartesian {
            x: -24934.457110145,
            y: -13883.157309732,
            z: -1662.925574181,
        };
        let b_cart = rtp.vec_to_cartesian(SphericalVec { coords: b_vec });
        assert_relative_eq!(b_cart.x, -25792.189012167, max_relative = 1.0e-8);
        assert_relative_eq!(b_cart.y, 9401.440618580, max_relative = 1.0e-8);
        assert_relative_eq!(b_cart.z, -7975.614708951, max_relative = 1.0e-8);
    }

    #[test]
    fn test_geod2geoc() {
        let geod = Geodetic {
            coords: Spherical {
                r: 1135.0,
                t: (90.0 - 45.5_f64).to_radians(),
                p: -23.5_f64.to_radians(),
            },
        };
        let geoc = geod.to_geocentric();
        assert_relative_eq!(geoc.coords.r, 1.177532882, max_relative = 1.0e-8);
        assert_relative_eq!(geoc.coords.t, 0.779521585, max_relative = 1.0e-8);
        assert_relative_eq!(geoc.coords.p, -0.410152374, max_relative = 1.0e-8);
    }

    #[test]
    fn test_geoc2geod() {
        let geod = Geocentric {
            coords: Spherical {
                r: 1.177532882,
                t: (90.0 - 45.336703245_f64).to_radians(),
                p: -23.500000111_f64.to_radians(),
            },
        };
        let geoc = geod.to_geodetic();
        assert_relative_eq!(geoc.coords.r, 1135.000000031, max_relative = 1.0e-8);
        assert_relative_eq!(
            geoc.coords.t,
            45.500000084_f64.to_radians(),
            max_relative = 1.0e-8
        );
        assert_relative_eq!(
            geoc.coords.p,
            -23.500000111_f64.to_radians(),
            max_relative = 1.0e-8
        );
    }
}
