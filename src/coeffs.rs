use crate::{AACGMv2Error, KMAX, NUM_COORDS, NUM_FLAGS, POLY_ORDER};
use chrono::{DateTime, Datelike, Timelike, Utc};
use std::path::Path;

/// Load a set of spherical harmonic coefficients.
///
/// Errors if the file is not a valid coefficient file (wrong number of coefficients or
/// coefficients cannot be interpreted as f64).
fn read_coeff_file<P: AsRef<Path>>(filename: P) -> Result<Vec<f64>, AACGMv2Error> {
    let coeffs_per_file = NUM_FLAGS * POLY_ORDER * NUM_COORDS * KMAX;

    // get the coefficients
    let contents = std::fs::read_to_string(filename)
        .map_err(|_| AACGMv2Error::CoeffFile("Could not read AACGM coefficients file"))?;

    let coeffs: Result<Vec<f64>, AACGMv2Error> = contents
        .split_whitespace()
        .map(|x| {
            x.parse::<f64>()
                .map_err(|_| AACGMv2Error::CoeffFile("Unable to interpret as f64"))
        })
        .collect();
    let flat_coeffs = coeffs?;
    if flat_coeffs.len() != coeffs_per_file {
        return Err(AACGMv2Error::CoeffFile(
            "Wrong number of coefficients in file",
        ));
    }

    Ok(flat_coeffs)
}

/// Load two sets of spherical harmonic coefficients.
///
/// Takes the 5-year epoch year prior to the desired time; bracketing set if +5 years.
///
/// Errors if `AACGM_v2_DAT_PREFIX` is not set in the environment, if the year is invalid,
/// or if the coefficient file is invalid.
fn load_coeffs(year: i32) -> Result<(Vec<f64>, Vec<f64>), AACGMv2Error> {
    /* default location of coefficient files */
    let prefix = match std::env::var("AACGM_v2_DAT_PREFIX") {
        Ok(x) => x,
        Err(_) => return Err(AACGMv2Error::Env("AACGM_v2_DAT_PREFIX")),
    };

    if year <= 0 {
        return Err(AACGMv2Error::CoeffFile("Invalid year for coefficients"));
    }

    let filename = format!("{prefix}{year}.asc");
    let coeffs = read_coeff_file(filename)?; /* forward coefficients */
    let filename = format!("{prefix}{}.asc", year + 5);
    let next_coeffs = read_coeff_file(filename)?; /* inverse coefficients */

    Ok((coeffs, next_coeffs))
}

/// Interpolate coefficients between adjacent 5-year epochs.
pub(crate) fn interpolate_coeffs(dt: DateTime<Utc>) -> Result<Vec<f64>, AACGMv2Error> {
    // epoch model year, round down to nearest 5 year increment
    let model_year = dt.year() / 5 * 5;

    /* fyear is the floating point time */
    let fdate = dt.ordinal0() as f64
        + (dt.hour() as f64 + (dt.minute() as f64 + dt.second() as f64 / 60.) / 60.) / 24.;
    let days_in_year = if dt.date_naive().leap_year() {
        366.0
    } else {
        365.0
    };
    let fyear = dt.year() as f64 + (fdate / days_in_year);

    let (forward, inverse) = load_coeffs(model_year)?;

    let mut sph_harm_model = vec![];
    /* time interpolation right here */
    for (x, y) in forward.iter().zip(inverse) {
        let val = x + (fyear - model_year as f64) * (y - x) / 5.0;
        sph_harm_model.push(val)
    }

    Ok(sph_harm_model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use chrono::NaiveDate;

    #[test]
    fn test_interpolate_coeffs() {
        let year = 2029;
        let month = 3;
        let day = 22;
        let hour = 3;
        let minute = 11;
        let second = 0;

        let dt = NaiveDate::from_ymd_opt(year, month, day)
            .unwrap()
            .and_hms_opt(hour, minute, second)
            .unwrap()
            .and_utc();
        let res = interpolate_coeffs(dt);
        assert!(res.is_ok());
        let x = res.unwrap();
        assert_relative_eq!(x[0], 0.11081712198135459);
    }
}
