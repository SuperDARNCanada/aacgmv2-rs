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

#[cfg(test)]
mod tests {
    use crate::aacgmv2::Aacgmv2;
    use crate::{Method, Transform, RADIUS_EARTH};
    use approx::assert_relative_eq;
    use chrono::NaiveDate;

    #[test]
    fn convert_coeffs() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let glat: f64 = 45.5;
        let glon: f64 = -23.5;
        let height: f64 = 1135.0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut aacgm_model = Aacgmv2::new(dt).unwrap();
        let cmd = Transform::GeodeticToAACGMv2;
        let flav = Method::Coeffs;
        let (out_lat, out_lon, out_rad) = aacgm_model
            .convert(glat, glon, height, &cmd, &flav)
            .unwrap();

        assert_relative_eq!(out_lat, 47.402897, max_relative = 1.0e-5);
        assert_relative_eq!(out_lon, 56.602300, max_relative = 1.0e-5);
        assert_relative_eq!(out_rad, 1.177533, max_relative = 1.0e-5);

        let cmd = Transform::AACGMv2ToGeodetic;
        let (glat, glon, height) = aacgm_model
            .convert(
                out_lat,
                out_lon,
                (out_rad - 1.0) * RADIUS_EARTH,
                &cmd,
                &flav,
            )
            .unwrap();

        assert_relative_eq!(glat, 45.439863, max_relative = 1.0e-5);
        assert_relative_eq!(glon, -23.477496, max_relative = 1.0e-5);
        assert_relative_eq!(height, 1134.977555, max_relative = 1.0e-5);
    }

    #[test]
    fn convert_trace() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let glat: f64 = 45.5;
        let glon: f64 = -23.5;
        let height: f64 = 1135.0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let mut aacgm_model = Aacgmv2::new(dt).unwrap();
        let cmd = Transform::GeodeticToAACGMv2;
        let flav = Method::Trace;
        let (out_lat, out_lon, out_rad) = aacgm_model
            .convert(glat, glon, height, &cmd, &flav)
            .unwrap();

        assert_relative_eq!(out_lat, 47.408678, max_relative = 1.0e-5);
        assert_relative_eq!(out_lon, 56.600154, max_relative = 1.0e-5);
        assert_relative_eq!(out_rad, 1.177533, max_relative = 1.0e-5);

        let cmd = Transform::AACGMv2ToGeodetic;
        let (glat, glon, height) = aacgm_model
            .convert(
                out_lat,
                out_lon,
                (out_rad - 1.0) * RADIUS_EARTH,
                &cmd,
                &flav,
            )
            .unwrap();

        assert_relative_eq!(glat, 45.500000, max_relative = 1.0e-5);
        assert_relative_eq!(glon, -23.500000, max_relative = 1.0e-5);
        assert_relative_eq!(height, 1135.000000, max_relative = 1.0e-5);
    }

    #[test]
    fn bad_date() {
        let year = 2030;
        let month = 1;
        let day = 1;
        let hour = 0;
        let minute = 0;
        let second = 0;

        let aacgm_model = Aacgmv2::new(
            NaiveDate::from_ymd_opt(year, month, day)
                .unwrap()
                .and_hms_opt(hour, minute, second)
                .unwrap()
                .and_utc(),
        );
        assert!(aacgm_model.is_err());
    }
}
