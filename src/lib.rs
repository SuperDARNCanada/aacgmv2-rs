#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

use thiserror::Error;

pub mod aacgmv2;
mod coeffs;
mod coords;
mod igrf;

const MAX_ALTITUDE: f64 = 2000.0; // km
const RADIUS_EARTH: f64 = 6371.2;
const NUM_FLAGS: usize = 2; // 0: geo to AACGM, 1: AACGM to geo
const NUM_COORDS: usize = 3; // x, y, z
const POLY_ORDER: usize = 5; // quartic polynomial fit in altitude
const SPH_HARM_ORDER: usize = 10; // order of spherical harmonic expansion
const KMAX: usize = (SPH_HARM_ORDER + 1) * (SPH_HARM_ORDER + 1); // number of spherical harmonic coefficients

#[derive(Error, Debug)]
pub enum AACGMv2Error {
    /// Invalid coordinates
    #[error("{0}")]
    Coords(String),

    /// Hit an error in an internal computation
    #[error("{0} - {1}")]
    Internal(i32, String),

    /// Hit an error in the IGRF module
    #[error("{0} - {1}")]
    Igrf(i32, String),

    /// Invalid environment
    #[error("{0}")]
    Env(&'static str),

    /// Invalid coefficients file
    #[error("{0}")]
    CoeffFile(&'static str),
}

/// How to conduct the coordinate transformation calculations.
pub enum Method {
    /// Use coefficients to calculate AACGM conversions
    Coeffs,

    /// Use field-line tracing to calculate AACGM conversions
    Trace,

    /// Use field-line tracing only above 2000 km
    AllowTrace,

    /// Use coefficients to calculate AACGM conversions, even above 2000km where coefficients are invalid
    BadIdea,
}

/// Coordinate transformation specifier.
#[derive(Clone, Debug)]
pub enum Transform {
    GeodeticToAACGMv2,
    AACGMv2ToGeodetic,
    GeocentricToAACGMv2,
    AACGMv2ToGeocentric,
}
