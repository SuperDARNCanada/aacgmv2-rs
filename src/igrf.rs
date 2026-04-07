use crate::coords::{Cartesian, Geographic, Magnetic, Spherical, SphericalVec};
use crate::{AACGMv2Error, RADIUS_EARTH};
use chrono::{DateTime, Datelike, Timelike, Utc};
use lazy_static::lazy_static;
use std::path::Path;

const IGRF_FIRST_EPOCH: i32 = 1590;
const IGRF_LAST_EPOCH: i32 = 2025;
const MAX_NUM_YEARS: usize = 100;
const IGRF_ORDER: usize = 13;
static IGRF_MAXK: usize = (IGRF_ORDER + 1) * (IGRF_ORDER + 1);

/// Struct containing the IGRF spherical harmonic coefficients, years, and secular variations.
struct IgrfCoeffSets {
    coeff_sets: [[f64; IGRF_MAXK]; MAX_NUM_YEARS],
    sec_vars: [f64; IGRF_MAXK],
}

#[derive(Clone, Debug, Default)]
pub(crate) struct IgrfCoeffs {
    coeffs: Vec<f64>,
    ecdip: Ecdip,
    geopack: Geopack,
    expansion_order: usize,
}

#[derive(Debug, Default, Copy, Clone)]
struct Ecdip {
    B02: f64,
    B0: f64,
    latref: f64,
    lonref: f64,
    L0: f64,
    L1: f64,
    L2: f64,
    E: f64,
    pos: [f64; 3],
    g2m: [[f64; 3]; 3],
}

#[derive(Debug, Default, Copy, Clone)]
struct Geopack {
    ctcl: f64,
    ctsl: f64,
    stcl: f64,
    stsl: f64,
    ct0: f64,
    st0: f64,
    cl0: f64,
    sl0: f64,
}

lazy_static! {
    static ref IGRF_COEFFS: IgrfCoeffSets = {
        // file containing the IGRF coefficients
    let filename = env!("IGRF_COEFFS");
        if filename.is_empty() {
        panic!("IGRF_COEFFS environment variable unset");
    }
        let (coeffs, sec_vars) = load_coeffs(filename).unwrap();
        IgrfCoeffSets { coeff_sets: coeffs, sec_vars }
    };
}

/// NAME:
///       IGRF_loadcoeffs
///
/// PURPOSE:
///       Load the entire set of spherical harmonic coefficients from the given
///       file.
///
///  Read the in the coefficients. Note that I am using the same ordering as
///  is used in the AACGM code. That is,
///
///   l    0  1  1  1  2  2  2  2  2  3  3  3  3  3  3  3  4  4  4  4  4 ...
///   m    0 -1  0  1 -2 -1  0  1  2 -3 -2 -1  0  1  2  3 -4 -3 -2 -1  0 ...
///
///  C & IDL index: k = l * (l+1) + m
///
///   k    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 ...
///
/// CALLING SEQUENCE:
///       err = IGRF_loadcoeffs();
///
/// Input Arguments:
///     filename      - name of file which contains IGRF coefficients; default
///                     is current IGRF model: igrf14coeffs.txt
///
/// Return Value:
///     Yearly and secular variation coefficients.
fn load_coeffs<P: AsRef<Path>>(
    filename: P,
) -> Result<([[f64; IGRF_MAXK]; MAX_NUM_YEARS], [f64; IGRF_MAXK]), AACGMv2Error> {
    let mut k: usize;
    let mut n: usize;
    let mut fac: f64;
    let mut iyear: usize;
    let mut nyear: usize;
    let mut dgrf = [0_i32; MAX_NUM_YEARS];
    let mut epoch: Vec<i32> = vec![];

    let mut sv: f64;
    let mut slm = [0.0_f64; IGRF_MAXK];
    let mut factorial = [0.0_f64; 2 * IGRF_ORDER + 1];
    let mut double_fac = [0.0_f64; 2 * IGRF_ORDER];

    let mut coeff_set = [[0.0_f64; IGRF_MAXK]; MAX_NUM_YEARS];
    let mut sec_vars = [0.0_f64; IGRF_MAXK];

    factorial[0] = 1.0;
    factorial[1] = 1.0;
    for k in 2..2 * IGRF_ORDER + 1 {
        factorial[k] = k as f64 * factorial[k - 1];
    }

    /* double factorial */
    double_fac[1] = 1.0;
    for k in 1..IGRF_ORDER {
        double_fac[2 * k + 1] = double_fac[2 * k - 1] * (2 * k + 1) as f64
    }

    for l in 0..IGRF_ORDER + 1 {
        for m in 0..l + 1 {
            k = l * (l + 1) + m; /* 1D index for (l, m) pair */
            n = l * (l + 1) - m; /* 1D index for (l, -m) pair */

            fac = if m != 0 { 2.0 } else { 1.0 };
            /* Davis 2004; Wertz 1978 recursion
            slm[k] = slm[n] = sqrt(fac*factorial[l-m]/factorial[l+m])*dfc[2*l-1]/factorial[l-m];
            */
            /* Winch 2004 */
            slm[k] = (fac * factorial[l - m] / factorial[l + m]).sqrt();
            slm[n] = slm[k]; // symmetric for +/- m
        }
    }

    // get the coefficients
    let contents = std::fs::read_to_string(filename)
        .map_err(|_| AACGMv2Error::Igrf(-1, "Could not read coefficients file".to_string()))?;

    let mut lines = contents.lines();

    // get third line
    let mut line = lines.nth(2).ok_or_else(|| {
        AACGMv2Error::Igrf(-1, "IGRF Coefficients file has too few lines".to_string())
    })?;

    nyear = 0;
    iyear = 0;
    for letter in line.chars() {
        match letter {
            'I' => {
                dgrf[iyear] = 0;
            }
            'D' => {
                dgrf[iyear] = 1;
            }
            'G' => {
                iyear += 1;
                nyear += 1;
            }
            _ => {}
        }
    }

    if nyear > MAX_NUM_YEARS {
        return Err(AACGMv2Error::Igrf(-2, "Too many years in file".to_string()));
    }

    /* get next line, which should have the following format:
     *
     * "g/h n m 1900.0 1905.0 ... 2010.0 2010-15"
     */
    line = lines.next().ok_or_else(|| {
        AACGMv2Error::Igrf(-1, "IGRF Coefficients file has too few lines".to_string())
    })?;

    // read the years, which should be 5-year integer epochs, skipping the secular variation column
    let mut years = line.split_whitespace();
    let _ = years.nth(2).ok_or_else(|| {
        AACGMv2Error::Igrf(
            -1,
            "IGRF Coefficients file year line improperly formatted".to_string(),
        )
    })?;

    for (m, year) in years.enumerate() {
        match year.parse::<f64>() {
            Ok(n) => {
                epoch.push(n.floor() as i32);
            }
            Err(_) => {
                if m < nyear {
                    return Err(AACGMv2Error::Igrf(
                        -2,
                        format!("Invalid year in file: {year}"),
                    ));
                }
            }
        }
    }

    // Now we are into the coefficients themselves
    for (i, line) in lines.enumerate() {
        let mut vals = line.split_whitespace();
        let sign_flag = vals
            .next()
            .ok_or_else(|| AACGMv2Error::Igrf(-1, "Empty coefficient line".to_string()))?;
        let sign: i32 = if sign_flag.contains('g') {
            1
        } else if sign_flag.contains('h') {
            -1
        } else {
            return Err(AACGMv2Error::Igrf(
                -1,
                format!("Invalid sign flag on coefficient line {i}"),
            ));
        };

        let degree = vals
            .next()
            .ok_or_else(|| {
                AACGMv2Error::Igrf(-1, format!("Coefficient line {i} missing `n` term"))
            })?
            .parse::<usize>()
            .map_err(|_| {
                AACGMv2Error::Igrf(-1, format!("Invalid degree `n` on coefficient line {i}"))
            })?;
        let order = vals
            .next()
            .ok_or_else(|| {
                AACGMv2Error::Igrf(-1, format!("Coefficient line {i} missing `m` term"))
            })?
            .parse::<usize>()
            .map_err(|_| {
                AACGMv2Error::Igrf(-1, format!("Invalid order `m` on coefficient line {i}"))
            })?;

        k = ((degree * (degree + 1)) as i32 + (order as i32 * sign)) as usize; // 1D index for (l, +/- m) pair,

        // Now iterate through coefficients for this (degree, order) combination by year
        let mut n: usize = 0;
        for val in vals {
            if n == nyear {
                // last column, secular variation
                sv = val.parse::<f64>().map_err(|_| {
                    AACGMv2Error::Igrf(
                        -1,
                        format!("Invalid secular variation for ({degree}, {order}) -> {val}"),
                    )
                })?;
                sec_vars[k] = sv * slm[k]; // Normalize
            } else if n < nyear {
                //
                let coeff = val.parse::<f64>().map_err(|_| {
                    AACGMv2Error::Igrf(
                        -1,
                        format!(
                            "Invalid coefficient for ({degree}, {order}), year {} -> {val}",
                            epoch[n]
                        ),
                    )
                })?;
                coeff_set[n][k] = coeff * slm[k]; // Normalize
            } else {
                return Err(AACGMv2Error::Igrf(
                    -1,
                    format!("Coefficient line {i} missing coefficients"),
                ));
            }
            n += 1;
        }
        if n != nyear + 1 {
            return Err(AACGMv2Error::Igrf(
                -1,
                "Not enough values on coefficient line".to_string(),
            ));
        }
    }

    Ok((coeff_set, sec_vars))
}

impl IgrfCoeffs {
    /// Interpolates coefficients for the given date.
    pub(crate) fn new(dt: DateTime<Utc>) -> Result<Self, AACGMv2Error> {
        let mut igrf_coeffs = IgrfCoeffs::default();
        let mut slm = [0.0_f64; IGRF_MAXK];
        let mut fctrl = [0.0_f64; 2 * IGRF_ORDER + 1];

        let mut coeffs = [0.0_f64; IGRF_MAXK];

        /* fyear is the floating point time */
        let fdate = dt.ordinal0() as f64
            + (dt.hour() as f64 + (dt.minute() as f64 + dt.second() as f64 / 60.) / 60.) / 24.;
        let days_in_year = if dt.date_naive().leap_year() {
            366.0
        } else {
            365.0
        };
        let fyear = dt.year() as f64 + (fdate / days_in_year);

        /* NOTE: FORTRAN code allows 10-year extrapolation beyond last epoch.
         * Here we are limiting to only 5 */
        if fyear < IGRF_FIRST_EPOCH as f64 || fyear > (IGRF_LAST_EPOCH + 5) as f64 {
            return Err(AACGMv2Error::Igrf(
                -3,
                "Cannot extrapolate to given datetime".to_string(),
            ));
        }

        let myear = dt.year() / 5 * 5; // epoch year, rounded down to nearest multiple of 5
        igrf_coeffs.expansion_order = if dt.year() < 1995 { 10 } else { 13 }; /* order of expansion */
        let i = (myear - IGRF_FIRST_EPOCH) as usize / 5; /* index of first set of coefs */

        let max_iter = (igrf_coeffs.expansion_order + 1) * (igrf_coeffs.expansion_order + 1);
        if fyear < IGRF_LAST_EPOCH as f64 {
            /* interpolate bounding coefficients */
            for (k, val) in coeffs.iter_mut().enumerate().take(max_iter).skip(1) {
                *val = IGRF_COEFFS.coeff_sets[i][k]
                    + (fyear - myear as f64)
                        * (IGRF_COEFFS.coeff_sets[i + 1][k] - IGRF_COEFFS.coeff_sets[i][k])
                        / 5.0;
            }
        } else {
            /* use secular variation */
            for (k, val) in coeffs.iter_mut().enumerate().take(max_iter).skip(1) {
                *val =
                    IGRF_COEFFS.coeff_sets[i][k] + (fyear - myear as f64) * IGRF_COEFFS.sec_vars[k];
            }
        }

        /* compute the components of the unit vector EzMag in geographic coordinates:
         * sin(theta0)*cos(lambda0), sin(theta0)*sin(lambda0)
         */

        /* C & IDL index: k = l * (l+1) + m */
        let g10 = -coeffs[2]; /* 1*2+0 = 2 */
        let g11 = coeffs[3]; /* 1*2+1 = 3 */
        let h11 = coeffs[1]; /* 1*2-1 = 1 */

        let sq = g11 * g11 + h11 * h11;

        let sqq = sq.sqrt();
        let sqr = (g10 * g10 + sq).sqrt();

        igrf_coeffs.geopack.sl0 = -h11 / sqq;
        igrf_coeffs.geopack.cl0 = -g11 / sqq;
        igrf_coeffs.geopack.st0 = sqq / sqr;
        igrf_coeffs.geopack.ct0 = g10 / sqr;

        igrf_coeffs.geopack.stcl = igrf_coeffs.geopack.st0 * igrf_coeffs.geopack.cl0;
        igrf_coeffs.geopack.stsl = igrf_coeffs.geopack.st0 * igrf_coeffs.geopack.sl0;
        igrf_coeffs.geopack.ctsl = igrf_coeffs.geopack.ct0 * igrf_coeffs.geopack.sl0;
        igrf_coeffs.geopack.ctcl = igrf_coeffs.geopack.ct0 * igrf_coeffs.geopack.cl0;

        /* for eccentric dipole coordinates */

        /* factorial for un-normalization */
        fctrl[0] = 1.0;
        fctrl[1] = 1.0;
        for k in 2..2 * IGRF_ORDER + 1 {
            fctrl[k] = k as f64 * fctrl[k - 1];
        }

        for l in 0..IGRF_ORDER + 1 {
            for m in 0..l + 1 {
                let k = l * (l + 1) + m; /* 1D index for l,m */
                let n = l * (l + 1) - m; /* 1D index for l,-m */

                let fac = if m != 0 { 2.0 } else { 1.0 };
                /* Davis 2004; Wertz 1978 recursion
                slm[k] = slm[n] = sqrt(fac*fctrl[l-m]/fctrl[l+m])*dfc[2*l-1]/fctrl[l-m];
                */
                /* Winch 2004 */
                slm[k] = (fac * fctrl[l - m] / fctrl[l + m]).sqrt();
                slm[n] = slm[k];
            }
        }

        /* S_(1,-1)^2 + S_(1,0)^2 + S_(1,1)^2 */
        igrf_coeffs.ecdip.B02 = coeffs[1] * coeffs[1] / (slm[1] * slm[1])
            + coeffs[2] * coeffs[2] / (slm[2] * slm[2])
            + coeffs[3] * coeffs[3] / (slm[3] * slm[3]);
        igrf_coeffs.ecdip.B0 = igrf_coeffs.ecdip.B02.sqrt();

        igrf_coeffs.ecdip.latref = (-coeffs[2] / slm[2] / igrf_coeffs.ecdip.B0)
            .asin()
            .to_degrees();
        igrf_coeffs.ecdip.lonref =
            180.0 + (coeffs[1] / slm[1]).atan2(coeffs[3] / slm[3]).to_degrees();

        let ca = igrf_coeffs.ecdip.latref.to_radians().cos();
        let sa = igrf_coeffs.ecdip.latref.to_radians().sin();
        let cb = igrf_coeffs.ecdip.lonref.to_radians().cos();
        let sb = igrf_coeffs.ecdip.lonref.to_radians().sin();

        igrf_coeffs.ecdip.g2m[0][0] = sa * cb;
        igrf_coeffs.ecdip.g2m[0][1] = sa * sb;
        igrf_coeffs.ecdip.g2m[0][2] = -ca;
        igrf_coeffs.ecdip.g2m[1][0] = -sb;
        igrf_coeffs.ecdip.g2m[1][1] = cb;
        igrf_coeffs.ecdip.g2m[1][2] = 0.;
        igrf_coeffs.ecdip.g2m[2][0] = cb * ca;
        igrf_coeffs.ecdip.g2m[2][1] = ca * sb;
        igrf_coeffs.ecdip.g2m[2][2] = sa;

        /*  2*S10*S20 + sqrt(3)*(S11*S21 + S1-1*S2-1)  */
        igrf_coeffs.ecdip.L0 = 2.0 * coeffs[2] / slm[2] * coeffs[6] / slm[6]
            + 3.0_f64.sqrt()
                * (coeffs[3] / slm[3] * coeffs[7] / slm[7]
                    + coeffs[1] / slm[1] * coeffs[5] / slm[5]);

        /* -S11*S20 + sqrt(3)*(S10*S21 + S11*S30 + S1-1*S2-2) */
        igrf_coeffs.ecdip.L1 = -coeffs[3] / slm[3] * coeffs[6] / slm[6]
            + 3.0_f64.sqrt()
                * (coeffs[2] / slm[2] * coeffs[7] / slm[7]
                    + coeffs[3] / slm[3] * coeffs[12] / slm[12]
                    + coeffs[1] / slm[1] * coeffs[4] / slm[4]);

        /* -S1-1*S20 + sqrt(3)*(S10*S2-1 - S1-1*S30 + S11*S2-2) */
        igrf_coeffs.ecdip.L2 = -coeffs[1] / slm[1] * coeffs[6] / slm[6]
            + 3.0_f64.sqrt()
                * (coeffs[2] / slm[2] * coeffs[5] / slm[5]
                    - coeffs[1] / slm[1] * coeffs[12] / slm[12]
                    + coeffs[3] / slm[3] * coeffs[4] / slm[4]);

        /* (L0*S10 + L1*S11 + L2*S1-1)/4/B02 */
        igrf_coeffs.ecdip.E = (igrf_coeffs.ecdip.L0 * coeffs[2] / slm[2]
            + igrf_coeffs.ecdip.L1 * coeffs[3] / slm[3]
            + igrf_coeffs.ecdip.L2 * coeffs[1] / slm[1])
            / 4.
            / igrf_coeffs.ecdip.B02;

        igrf_coeffs.ecdip.pos[0] = RADIUS_EARTH
            * (igrf_coeffs.ecdip.L1 - coeffs[3] / slm[3] * igrf_coeffs.ecdip.E)
            / 3.
            / igrf_coeffs.ecdip.B02;
        igrf_coeffs.ecdip.pos[1] = RADIUS_EARTH
            * (igrf_coeffs.ecdip.L2 - coeffs[1] / slm[1] * igrf_coeffs.ecdip.E)
            / 3.
            / igrf_coeffs.ecdip.B02;
        igrf_coeffs.ecdip.pos[2] = RADIUS_EARTH
            * (igrf_coeffs.ecdip.L0 - coeffs[2] / slm[2] * igrf_coeffs.ecdip.E)
            / 3.
            / igrf_coeffs.ecdip.B02;

        igrf_coeffs.coeffs = coeffs.to_vec();

        Ok(igrf_coeffs)
    }

    /// Calculate the magnetic field at a point.
    pub(crate) fn compute(&self, point: &Spherical) -> Result<SphericalVec, AACGMv2Error> {
        let mut tbrtp = [0.0_f64; 3];
        let mut brtp = Cartesian::default();
        let mut cosm_arr = [0.0_f64; IGRF_ORDER + 1];
        let mut sinm_arr = [0.0_f64; IGRF_ORDER + 1];

        // Must avoid singularity at the poles (dividing by sin(theta) later)
        let st = point.t.sin();
        let dt = if st < 0.0 { 1.0e-15 } else { -1.0e-15 };
        let theta = if st.abs() < 1e-15 {
            point.t + dt
        } else {
            point.t
        };

        /* Compute the values of the Legendre Polynomials, and derivatives */
        let (plm_val, plm_val_deriv) = igrf_plm(theta, self.expansion_order)?;

        let aor = 1. / point.r; // r is in units of RE to be consistent with geopack
        let mut afac = aor * aor;

        // array of trig functions in phi for faster computation
        for k in 0..IGRF_ORDER + 1 {
            cosm_arr[k] = (k as f64 * point.p).cos();
            sinm_arr[k] = (k as f64 * point.p).sin();
        }

        for l in 1..self.expansion_order + 1 {
            // no l = 0 term in IGRF
            tbrtp[0] = 0.0;
            tbrtp[1] = 0.0;
            tbrtp[2] = 0.0;
            for m in 0..l + 1 {
                let k = l * (l + 1) + m; // g
                let n = l * (l + 1) - m; // h

                tbrtp[0] +=
                    (self.coeffs[k] * cosm_arr[m] + self.coeffs[n] * sinm_arr[m]) * plm_val[k];
                tbrtp[1] += (self.coeffs[k] * cosm_arr[m] + self.coeffs[n] * sinm_arr[m])
                    * plm_val_deriv[k];
                tbrtp[2] += (-self.coeffs[k] * sinm_arr[m] + self.coeffs[n] * cosm_arr[m])
                    * m as f64
                    * plm_val[k];
            }
            afac *= aor;

            brtp.x += afac * (l + 1) as f64 * tbrtp[0];
            brtp.y -= afac * tbrtp[1];
            brtp.z -= afac * tbrtp[2];
        }

        brtp.z /= theta.sin();

        Ok(SphericalVec { coords: brtp })
    }

    /// Convert geographic coordinates to magnetic.
    pub(crate) fn geo2mag(&self, geo: Geographic) -> Magnetic {
        let mut mag = Magnetic::default();
        mag.coords.x = geo.coords.x * self.geopack.ctcl + geo.coords.y * self.geopack.ctsl
            - geo.coords.z * self.geopack.st0;
        mag.coords.y = geo.coords.y * self.geopack.cl0 - geo.coords.x * self.geopack.sl0;
        mag.coords.z = geo.coords.x * self.geopack.stcl
            + geo.coords.y * self.geopack.stsl
            + geo.coords.z * self.geopack.ct0;

        mag
    }

    /// Convert magnetic coordinates to geographic.
    pub(crate) fn mag2geo(&self, mag: Magnetic) -> Geographic {
        let mut geo = Geographic::default();

        geo.coords.x = mag.coords.x * self.geopack.ctcl - mag.coords.y * self.geopack.sl0
            + mag.coords.z * self.geopack.stcl;
        geo.coords.y = mag.coords.x * self.geopack.ctsl
            + mag.coords.y * self.geopack.cl0
            + mag.coords.z * self.geopack.stsl;
        geo.coords.z = mag.coords.z * self.geopack.ct0 - mag.coords.x * self.geopack.st0;

        geo
    }
}

/*-----------------------------------------------------------------------------
;
; NAME:
;       IGRF_Plm
;
; PURPOSE:
;       Internal function to compute array of Gaussian Normalized Associated
;       Legendre functions and the corresponding derivatives.
;
; CALLING SEQUENCE:
;       err = IGRF_Plm(theta, order, plmval, dplmval);
;
;     Input Arguments:
;       theta         - co-latitude in radians
;       order         - order of expansion, should NOT exceed IGRF_ORDER
;
;     Output Arguments:
;       plmval        - pointer to array for storage of values
;       dplmval       - pointer to array for storage of derivative values
;
;     Return Value:
;       error code
;
;     Notes: I am using array indexing similar to that used for m=-l to l,
;            but here m=0 to l, so the arrays are too big and there are no
;            values stored in locations for m<0. Probably should fix that...
;
;       values are stored in a 1D array of dimension (order+1)^2. The
;       indexing scheme used is:
;
;             g  h  g  g  h  h  g  g  g  h  h  h  g  g  g  g  h  h  h  h  h ...
;        l    0  1  1  1  2  2  2  2  2  3  3  3  3  3  3  3  4  4  4  4  4 ...
;        m    0 -1  0  1 -2 -1  0  1  2 -3 -2 -1  0  1  2  3 -4 -3 -2 -1  0 ...
;C & IDL j    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 ...
;FORTRAN j    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 ...
;
;+-----------------------------------------------------------------------------
*/
fn igrf_plm(
    theta: f64,
    order: usize,
) -> Result<([f64; IGRF_MAXK], [f64; IGRF_MAXK]), AACGMv2Error> {
    let mut plmval = [0.0_f64; IGRF_MAXK];
    let mut dplmval = [0.0_f64; IGRF_MAXK];

    if order > IGRF_ORDER {
        return Err(AACGMv2Error::Igrf(-1, format!("Cannot compute Legendre polynomials with order above IGRF_ORDER ({order} > {IGRF_ORDER})")));
    };

    let st = theta.sin();
    let ct = theta.cos();

    plmval[0] = 1.; // l=0, m=0
    dplmval[0] = 0.; // l=0, m=0

    // compute values of P^{l,l} and dP^{l,l}/dtheta
    for l in 1..order + 1 {
        let k = l * (l + 1) + l; // l = m
        let n = (l - 1) * l + l - 1; // l-1 = m-l, i.e., previous l=m
        let a = (2 * l - 1) as f64;
        plmval[k] = a * plmval[n] * st;
        dplmval[k] = a * (dplmval[n] * st + plmval[n] * ct);
    }

    plmval[2] = ct; // 1,0
    dplmval[2] = -st; // 1,0
                      // compute values of P^{l,m} and dP^{l,m}/dtheta
    for l in 2..order + 1 {
        for m in 0..l {
            let k = l * (l + 1) + m; /* l,m */
            let n = (l - 1) * l + m; /* l-1,m */
            let p = (l - 2) * (l - 1) + m; /* l-2,m */

            // numerical recipes in C
            let a = (2 * l - 1) as f64;
            if m == l - 1 {
                plmval[k] = a * ct * plmval[n] / (l - m) as f64;
                dplmval[k] = a * (ct * dplmval[n] - st * plmval[n]) / (l - m) as f64;
            } else {
                let b = (l + m - 1) as f64;
                plmval[k] = (a * ct * plmval[n] - b * plmval[p]) / (l - m) as f64;
                dplmval[k] =
                    (a * (ct * dplmval[n] - st * plmval[n]) - b * dplmval[p]) / (l - m) as f64;
            }
        }
    }

    Ok((plmval, dplmval))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use chrono::NaiveDate;

    #[test]
    fn test_igrf() {
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
        let model = IgrfCoeffs::new(dt).unwrap();

        let point = Spherical {
            r: 1.177532887,
            t: 0.779521583,
            p: 5.873032931,
        };
        let res = model.compute(&point);
        assert!(res.is_ok());
        let rtp = res.unwrap();
        assert_relative_eq!(rtp.radial(), -24934.457110145, max_relative = 1.0e-8);
        assert_relative_eq!(rtp.colatitude(), -13883.157309732, max_relative = 1.0e-8);
        assert_relative_eq!(rtp.longitude(), -1662.925574181, max_relative = 1.0e-8);
    }

    #[test]
    fn test_mag2geo() {
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
        let model = IgrfCoeffs::new(dt).unwrap();

        let mag = Magnetic {
            coords: Cartesian {
                x: 1.201891469,
                y: 1.822774737,
                z: 0.0,
            },
        };
        let geo = model.mag2geo(mag);
        assert_relative_eq!(geo.coords.x, 2.090921575, max_relative = 1.0e-8);
        assert_relative_eq!(geo.coords.y, -0.599541509, max_relative = 1.0e-8);
        assert_relative_eq!(geo.coords.z, -0.188806233, max_relative = 1.0e-8);
    }
}
